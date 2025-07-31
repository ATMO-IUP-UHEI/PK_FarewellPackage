# skript to select ObsPack data for regional Flexpart Inversion of North American CO2 fluxes

# select region, timeperiod
# only use in-situ and tower measurements, for towers with multiple inlet height, only use highest inlet
# only use with CarbonTracker assimilation flag CT_assim=1


import numpy as np
import xarray as xr
import pandas as pd
import os
import datetime as dt
import re
from collections import defaultdict
from timezonefinder import TimezoneFinder
import pytz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

# Initialize timezone finder
tf = TimezoneFinder()

def convert_to_local_time(utc_time, lat, lon):
    """Convert UTC to local time based on latitude and longitude.
    Args:
        utc_time (dt.datetime): UTC time
        lat (float): Latitude.
        lon (float): Longitude.
    Returns:
        UTC datetime and timezone_str
        or None if timezone not found
    """
    timezone_str = tf.timezone_at(lng=lon, lat=lat)
    if timezone_str:
        local_tz = pytz.timezone(timezone_str)
        return (utc_time.tz_convert(local_tz).replace(tzinfo=None), timezone_str)
    return None  # Return None if timezone is not found

def convert_to_utc(local_time, lat, lon):
    """Convert local date and time to UTC based on latitude and longitude.
    Args:
        local_time (dt.datetime): Local time as datetime
        lat (float): Latitude.
        lon (float): Longitude.
    Returns:
        UTC datetime or None if timezone not found.
    """
    timezone_str = tf.timezone_at(lng=lon, lat=lat)
    if timezone_str:
        local_tz = pytz.timezone(timezone_str)
        # Combine date and time into a datetime object
        local_dt = local_tz.localize(local_time)
        # Convert to UTC
        utc_dt = local_dt.astimezone(pytz.utc)
        return utc_dt.replace(tzinfo=None)
    return None  # Timezone not found

def sel_stations(path,region, start_date, end_date, s_dir):
    '''  Select surface and tower stations in selected region and time period, only use tower files with max intake height
    Args:
        path: path to folder containing the obspack data
        region: list [lat_min, lat_max, lon_min, lon_max] for region that should be selected
        start_date, end_date: dt.date objects giving the time period
        s_dir: path to directory where selected files should be stored
    Returns: 
        nothing, saves selected station files into s_dir
    '''
    # list of surface and tower files
    files=[f for f in os.listdir(path) if ('surface' in f or 'tower' in f)]
    files.sort()

    # only use tower files with max intake height
    # Dictionary to store the max height file for each prefix
    file_groups = defaultdict(lambda: ("", 0))
    other_files=[]
    for file in files:
        # print(file)
        match = re.match(r"(.+)-(\d+)magl\.nc", file)
        if match:
            prefix, height = match.groups()
            height = int(height)
            if height > file_groups[prefix][1]:  # Keep max height
                file_groups[prefix] = (file, height)
        else: # if no matching height given, keep file
            other_files.append(file)

    # Get only the files with max height
    filtered_files = [file for file, _ in file_groups.values()] + other_files
    filtered_files.sort()
    
    # for all filtered files, save with local time
    # get desired time range and region for filtered files
    for f in filtered_files:
        ds=xr.open_dataset(path+f)
        # ds=ds.assign_coords(time=ds.time, latitude=ds.latitude, longitude=ds.longitude, altitude=ds.altitude)
        temp=ds.drop_dims(['calendar_components','dim_concerns' ])
        # check that there are measurements with CT_assim=1 flag for desired timeperiod & region with 
        if temp.where((temp.CT_assim==1) & (temp.time.dt.date>=start_date)&(temp.time.dt.date<=end_date) & (temp.latitude> region[0])&(temp.latitude<region[1]) & (temp.longitude> region[2])&(temp.longitude<region[3]), drop=True).sizes['obs']>0:
            ds_sel=ds.where((ds.CT_assim==1) & (ds.time.dt.date>=start_date)&(ds.time.dt.date<=end_date) & (ds.latitude> region[0])&(ds.latitude<region[1]) & (ds.longitude> region[2])&(ds.longitude<region[3]), drop=True)
            # get local time and timezone
            local_times, timezones = zip(*np.array([convert_to_local_time(t, lat, lon) for t, lat, lon in zip(pd.to_datetime(ds_sel.time.values, utc=True), ds_sel.latitude.values, ds_sel.longitude.values) ]))
            ds_sel['local_time']=(('obs'), np.array(local_times, dtype='datetime64[ns]'))
            ds_sel['local_time'].attrs = {
                'description': 'Local time converted from UTC based on lat/lon. Time-averaged values are reported at the beginning of the averaging interval.',
                'source': 'TimezoneFinder and pytz'
            }
            ds_sel['timezone']=(('obs'), np.array(timezones))
            # save dile as .nc file
            ds_sel.to_netcdf(s_dir+f)
            print(f"saved to {s_dir}{f}")
    return

def get_4h_mean_stations(dir_path, sdir_path, h_lim=1100):
    ''' Get daily 4h mean values for each station (local afternoon 4h mean for surface meas, local night 4h mean for mountaintop stations)
    Args:
        dir_path: path to directory containing selected data
        sdir_path: path to directory where data should be saved
        hlim: elevation above which station is considered a mountain station, defaults to 1100masl (m above sea level)
    Returns:
        nothing, saves daily 4h mean values for each station file
    '''
    files=[f for f in os.listdir(dir_path) if f.endswith('.nc')]
    for f in files:
        data=xr.open_dataset(dir_path+f)
        if len(np.unique(data.elevation.values))!=1:
            print(f'Error with elevation for file {f}, should be exactly one value!')
            print(data.elevation.values)
        else:
            elevation=np.unique(data.elevation.values).item()
            if elevation<h_lim: # surface meas, local afternoon 4hour mean
                start_time=dt.time(12,0)
                end_time=dt.time(16,0)          
            elif elevation>h_lim:   # mountain, nightly 4hour mean
                start_time=dt.time(0,0)
                end_time=dt.time(4,0)
            data_sel=data.where((data.local_time.dt.time>=start_time)&(data.local_time.dt.time<end_time))
            # if not all values are nans
            if not np.isnan(data_sel.value).all().item():
                # drop nans
                data_sel=data.where((data.local_time.dt.time>=start_time)&(data.local_time.dt.time<end_time), drop=True)    # do it like this instead of using drop nan because of oter nans in data, would drop to much
                if not 'nvalue' in data_sel.keys():    # if number of meas values per average not given, assume its one
                    data_sel['nvalue']=(('obs'), np.ones(data_sel.sizes['obs']))
                data_mean=data_sel.groupby("local_time.date").mean(dim='obs')    # mean for local 4 hour period
                data_mean['num_meas']=data_sel.groupby("local_time.date").count(dim='obs').nvalue  
                data_mean['num_meas']=data_mean['num_meas'].assign_attrs(description='Number of hourly measurements used for 4h mean value')    # can be larger than 4 if hourly averages given e.g. every 30 minutes ('co2_bao_tower-insitu_1_allvalid-300magl.nc')
                data_mean['num_toal']=data_sel.groupby("local_time.date").sum(dim='obs').nvalue
                data_mean['num_toal']=data_mean['num_toal'].assign_attrs(description='Sum over Number of individual measurements used for 4h mean value')
                # local time at start of 4h averaging period
                data_mean['local_time']=(('date'),[dt.datetime.combine(d, start_time) for d in data_mean.date.values])
                # get start_time in UTC
                data_mean['time'] = (('date'), [convert_to_utc(t, lat, lon) for t, lat, lon in zip(pd.to_datetime(data_mean.local_time.values), data_mean.latitude.values, data_mean.longitude.values)])
                data_mean['time']=data_mean['time'].assign_attrs(description='UTC values are reported at the beginning of the 4hour averaging interval.')
                data_mean=data_mean.swap_dims({"date": "time"}).drop_vars('date')
                data_mean.to_netcdf(sdir_path+f)
                print(f'saved to {sdir_path+f}')
    return

def plot_meas_station_map(path, region, sfig_path='is_meas_map.png'):
    ''' Plot figure with in-situ measurement station locations within selected region.
    Args:
        path: path to Obspack data with one file per station
        region: [lat_min, lat_max, lon_min, lon_max]
        sfig_path: path to where figure should be saved, defaults to 'is_meas_map.png'
    Returns: nothig, saves figure'''
    
    files=[f for f in os.listdir(path)]
    # distinguish between surface and tower, insitu and flask
    surf_is=[f for f in files if 'surface-insitu' in f]
    surf_pfp=[f for f in files if 'surface-pfp' in f]
    surf_flask=[f for f in files if 'surface-flask' in f]   # all flask meas
    tower_is=[f for f in files if 'tower-insitu' in f]

    fig, ax=plt.subplots(1,1,figsize=(12,8), subplot_kw={'projection': ccrs.PlateCarree()})
    for f in surf_is:
        temp=xr.open_dataset(path+f)
        ax.scatter(temp.longitude[0], temp.latitude[0], marker='o', color='royalblue')
    for f in surf_flask:
        temp=xr.open_dataset(path+f)
        ax.scatter(temp.longitude[0], temp.latitude[0], marker='^', color='darkorange')
    for f in surf_pfp:
        temp=xr.open_dataset(path+f)
        ax.scatter(temp.longitude[0], temp.latitude[0], marker='d', color='#f781bf')
    for f in tower_is:
        temp=xr.open_dataset(path+f)
        ax.scatter(temp.longitude[0], temp.latitude[0], marker='s', color='#4daf4a')

    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    # Add coastlines
    ax.coastlines(linewidth=0.5, alpha=0.3)
    # Add U.S. state borders
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5,alpha=0.3)
    gl = ax.gridlines(draw_labels=True,color='gray', alpha=0.3, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 20}
    gl.ylabel_style = {'size': 20}
    ax.grid(c='lightgrey')
    ax.set_xlim(region[2],region[3]) 
    ax.set_ylim(region[0],region[1])
    ax.tick_params(labelsize=30)  # Increase lat/lon label size

    # Manually create legend handles
    custom_legend = [
        Line2D([0], [0], marker='o', label='surface insitu',color='royalblue', linestyle=''),
        Line2D([0], [0], marker='^', label='surface flask',color='darkorange', linestyle=''),
        Line2D([0], [0], marker='d', label='surface pfp',color='#f781bf', linestyle=''),
        Line2D([0], [0], marker='s', label='tower insitu',color='#4daf4a', linestyle='')
    ]
    # Add custom legend
    ax.legend(handles=custom_legend, loc='lower left', fontsize=20)
    # fig.tight_layout()
    plt.show()
    fig.savefig(path+sfig_path, dpi=300, bbox_inches='tight')
    return

#########################################################################################

if __name__ == '__main__':
    flexpart_region=[12,56,-134,-62]  # lat_min, lat_max, lon_min, lon_max
    start_date=dt.date(2010,6,1)
    end_date=dt.date(2010,6,30)      # including end_date
    # path to obspack data
    obspack_path='/net/dsvr-02/mnt/data2/users/eschoema/ObsPack/obspack_co2_1_GLOBALVIEWplus_v5.0_2019-08-12/data/nc/'
    
    # select time period, region, highest towe inlet
    sel_sdir='/home/pkuehn/data/ObsPack/test/test_sel/'
    if not os.path.isdir(sel_sdir):
        os.makedirs(sel_sdir)
    sel_stations(obspack_path,flexpart_region, start_date, end_date, sel_sdir)
    # plot all selected ObsPack measurement stations
    plot_meas_station_map(sel_sdir, flexpart_region)
    
    # calculate 4h mean
    # mountain for elevation > 1100 masl
    h_lim=1100  #masl
    # path to save 4h mean values
    mean_sdir='/home/pkuehn/data/ObsPack/test/test_mean/'
    if not os.path.isdir(mean_sdir):
        os.makedirs(mean_sdir)
    get_4h_mean_stations(sel_sdir, mean_sdir, h_lim)
    
    # combine individual station files to one dataset
    spath=f'/home/pkuehn/data/ObsPack/test/test_lat{flexpart_region[0]}_{flexpart_region[1]}_lon{flexpart_region[2]}_{flexpart_region[3]}sel_mean_combined.nc'
    files=[f for f in os.listdir(mean_sdir) if f.endswith('.nc')]
    ds_list=[]
    for f in files:
        ds=xr.open_dataset(mean_sdir+f).expand_dims(file=[f])
        ds["num_meas"] = ds["num_meas"].astype(float)
        ds["local_time"].encoding["dtype"] = "float32"  # Use a floating-point type
        ds["local_time"].encoding["units"] = "hours since 2009-10-06T12:00:00"
        ds_list.append(ds)
    ds=xr.concat(ds_list, dim='file')
    ds.to_netcdf(spath)
    return


import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import os
import re
from dateutil.relativedelta import relativedelta

# get monday of the week containing the date
def get_start_date_of_week(date):
    start_date=date-dt.timedelta(days=date.weekday())   # dt.date().weekday() goes from 0 to 6
    return start_date

def get_scaling_factors(flexpart_region, start_date, end_date, dir_path, spath):
    ''' Caclulate diurnal cycle scaling factors based on 3hourly TM5-4DVar prior data
    Args:
        dir_path: path to data directory with 3hrly prior data
        flexpart_region: list of lat_min, lat_max, lon_min, lon_max that the dataset should be cut to
        start_date, end_date: dt.date objects defining the time period of flexpart releases
        spath: path to where data should be saved
    Retruns:
        nothing, saves data as .nc file
    '''
    # preprocessing, cut region
    lat_list=np.arange(flexpart_region[0]+0.5, flexpart_region[1], step=1)
    lon_list=np.arange(flexpart_region[2]+0.5, flexpart_region[3], step=1)
    # get start of flexpart simulations, duration is 10 days, so 10 days before earliest release
    sim_start=get_start_date_of_week(start_date+relativedelta(days=-10))

    def preprocess(ds):
        # preprocess function that reads the date from filename, gets datetime from timestep and 
        # cuts to flexpart_region
        # calculates scaling factor
        
        # Extract the filename from the internal attribute
        filename = os.path.basename(ds.encoding['source'])  # get filename from dataset encoding
        match = re.search(r'(\d{8})', filename)  # extract date part (YYYYMMDD)
        if match:
            date_str = match.group(1)
            date = pd.to_datetime(date_str, format='%Y%m%d')
            ds = ds.expand_dims({'date': [date]})  # add 'date' dimension
        else:
            raise ValueError(f"Could not find date in filename {filename}")
        # drop unnecessary variables
        ds=ds.drop_vars(['nee_total', 'fire_total', 'fossil_total', 'ocean_total'])
        # get latitude, longitude
        ds=ds.assign_coords(latitude=(('nlat'), ds.lat_edges.values[:-1]+1), longitude=(('nlon'), ds.lon_edges.values[:-1]+1.5))
        # Cut to flexpart region
        ds=ds.swap_dims({'nlat':'latitude', 'nlon':'longitude'})
        # regid onto 1x1 grid
        ds=ds.interp(latitude=lat_list,longitude=lon_list,method='nearest').sel(latitude=slice(flexpart_region[0], flexpart_region[1]), longitude=slice(flexpart_region[2], flexpart_region[3]))
        # get time in UTC
        ds=ds.assign_coords(time=(('time_step'), (ds.time_step.values*3+1.5), {'description':'UTC time of day [h]'}))   #date+pd.to_timedelta((ds.time_step.values*3+1.5)
        ds=ds.swap_dims({'time_step':'time'}).drop_vars(['lat_edges', 'lon_edges'])
        return ds

    # Build the list of filepaths
    file_list = [f'{dir_path}flux1x1_{date.strftime("%Y%m%d")}.nc' for date in pd.date_range(sim_start, end_date)]
    # sort file list
    file_list.sort()
    
    # combine into one ds
    # Now use open_mfdataset
    data = xr.open_mfdataset(
        file_list,
        preprocess=preprocess,
        combine='by_coords',
        chunks="auto",  # Open lazily with dask
        parallel=True
    )
    
    # get weekly mean
    date_min=pd.to_datetime(data.date.min().values)
    date_max=pd.to_datetime(data.date.max().values)
    # Create the 7-day period bins
    period_bins = pd.date_range(start=get_start_date_of_week(date_min), end=(get_start_date_of_week(date_max)+dt.timedelta(days=7)), freq="7D")
    # Cut the time array into the defined 7-day bins
    time_bins = pd.cut(data["date"], bins=period_bins, right=False, labels=period_bins[:-1])
    # Assign the new time bins as coordinates
    data=data.assign_coords(date=time_bins)
    # weekly mean
    # the timestamp corresponds to the first day of the week
    data=data.groupby("date").mean('date')
    for col in ['fossil', 'nee', 'ocean', 'fire']:
        data[col]=data[col].assign_attrs(description='mean weekly flux')
    # get scaling of biospheric flux 
    data['nee_mean']=data.nee.mean(dim='time').assign_attrs(description='daily mean of weekly flux', units='gC/m^2/day')
    # if mno nee flux, would divide by zero, so set scaling factor to zero 
    data['scaling_bio'] = xr.where(data.nee_mean != 0,data.nee / data.nee_mean,0).assign_attrs(description='scaling of nee fluxes, nee/nee_weekly_mean', units='[]')
    # data['scaling_bio']=(data.nee/data.nee_mean).assign_attrs(description='scalling of nee fluxes, nee/nee_weekly_mean', units='[]')
    # get total flux
    data['total_flux']=(data.nee+data.fire+data.fossil+data.ocean).assign_attrs(description='total mean weekly flux', units='gC/m^2/day')
    data['total_flux_mean']=data.total_flux.mean(dim='time').assign_attrs(description='daily mean of total weekly flux', units='gC/m^2/day')
    # get rest flux, everything except nee
    data['rest_flux']=(data.fire+data.fossil+data.ocean).assign_attrs(description='mean of rest weekly flux, fossil+fire+ocean', units='gC/m^2/day')
    data['rest_flux_mean']=data.rest_flux.mean(dim='time').assign_attrs(description='daily mean of rest weekly flux', units='gC/m^2/day')
    # get total scaling factor
    # sc_tot = sc_bio*bio_mean/tot_mean + rest_mean/tot_mean
    # sc_tot = bio_h/tot_mean + rest_mean/tot_mean      with sc_bio = bio_h/bio_mean
    data['scaling_total']=(data.nee/data.total_flux_mean+data.rest_flux_mean/data.total_flux_mean).assign_attrs(description='scalling of total fluxes, sc_tot = bio_h/tot_mean + rest_mean/tot_mean', units='[]')

    # Convert weekly data to daily by forward-filling values
    data = data.reindex(date=data.date, method='ffill')
    # save data
    data.to_netcdf(spath)
    print(f'saved to {spath}')
    return
def interpolate_total_scaling(path, spath):
    ''' Interpolate total scaling factors from mean weekly to hourly values
    Args: 
        path: path to mean weekly TM5-4DVar prior scalings
        spath: path to where output data should be saved
    '''
    data=xr.open_dataset(path)
    # Shift time=1.5 to time=25.5 (1.5 + 24) to close the circl and combine
    data_interp=xr.concat([data, data.sel(time=1.5).assign_coords(time=25.5)], dim='time')
    # 4. Create new finer time grid (hourly, 0.5 to 23.5, but allow for 24+ values)
    new_time = np.arange(1.5, 25.5, 1.0)  # 1-hour steps
    # 5. Interpolate
    data_interp = data_interp[['scaling_total']].interp(time=new_time, method='linear')
    # 6. Modulo time back into [0, 24)
    data_interp = data_interp.assign_coords(time=(data_interp.time % 24)).sortby('time')

    # forward fill dates
    # Convert monthly data to daily by forward-filling values
    start_date=data_interp.date.min().dt.date.item()
    end_date=data_interp.date.max().dt.date.item()+relativedelta(days=6)
    data_interp = data_interp.reindex(date=pd.date_range(start_date,end_date), method='ffill')

    # get times as datetimes, need same times as footprint data
    # 1. Convert hours to timedelta
    hour_offset = xr.apply_ufunc( pd.to_timedelta, data_interp['time']-0.5, kwargs={'unit': 'h'}, vectorize=True)
    # 2. Add to date
    datetime = data_interp.date + hour_offset

    # stack date and time dimensions
    data_interp=data_interp.stack(datetime_index=('date', 'time')).drop_vars(['date', 'time'])
    # add datetime with name 'time', so its the same as for the footprints
    data_interp['time']=(('datetime_index'), datetime.values.flatten())
    data_interp=data_interp.swap_dims({'datetime_index':'time'})
    data_interp.to_netcdf(spath)
    return

def CutAndGetPressuresTM5mix(filepath, region, spath, daystr):
    [Lat_min, Lat_max, Long_min,Long_max]=region
    # read presssure params a & b
    p_param=xr.open_dataset(filepath)
    # read co2 data & surface pressure
    data=xr.open_dataset(filepath, group = 'glb300x200')
    # select box region
    data=data.sel(latitude=slice(Lat_min, Lat_max), longitude=slice(Long_min,Long_max)).drop_vars(['lat_edges','lon_edges'])  # drop unnecessary labels
    # get times as datetime
    data['times']=np.array([dt.datetime(*row) for row in data.sample_times.values])

    # calculate pressure at boundaries
    data['p_boundary']=p_param.at+p_param.bt*data.pressure
    # get pressure difference for each level
    data['p_diff']=(data.p_boundary.sel(boundaries=slice(0,25))-data.p_boundary.sel(boundaries=slice(1,26))).rename({'boundaries': 'levels'})

    # get total column xco2, pressure weighed using p_diff
    # for each level weigh mixing ratio with pressure difference over that level, sum over levels, divide by pressure difference over entire column
    data['xco2']=((data.mix*data.p_diff).sum(dim='levels'))/(data.p_boundary.sel(boundaries=0)-data.p_boundary.sel(boundaries=25))
    data['xco2_mean']=data.xco2.mean(dim='times')
    
    # save dataset 
    data.to_netcdf(spath+f'xco2_mean_{daystr}.nc')
    return

# main
if __name__ == "__main__": 
    # get diurnal cycle scaling factors
    flexpart_region=[12,56,-134,-62]       # lat_min, lat_max, lon_min, lon_max
    # read data for entire period
    start_date=dt.date(2009,10,1)
    end_date=dt.date(2011,3,31)
    # path to 3hrly prior data
    dir_path='/mnt/data/users/eschoema/TM5Inversion/3h_fluxes/three_hourly_flux/RemoTeC+IS (lo)/apri/'
    # path to save 3hrly weakly mean scaling
    scaling_path=f'/home/pkuehn/data/TM5Inversion/weekly_mean_prior_scaling_RemoTeC+IS.nc'
    get_scaling_factors(flexpart_region, start_date, end_date, dir_path, scaling_path)
    # path to save hourly interpolated scaling factors
    spath="/home/pkuehn/data/TM5Inversion/high_res_total_scaling_RemoTeC+IS.nc"
    interpolate_total_scaling(scaling_path, spath)
    
    # get cut molefraction fields
    # TODO check which are bias corected
    Dtypes=['RemoTeC+IS', 'ACOS+IS', 'IS']  # 'RemoTeC+IS', 'ACOS+IS', 'IS', 'RemoTeC+IS-land_ocean_bc'
    # (old) flexpart_region=[6,62,-140,-56]       # lat_min, lat_max, lon_min, lon_max
    # box region, extend two grid cells around, 2x3° grid
    Lat_min,Lat_max=2,66
    Long_min,Long_max=-146, -50
    # select year and month
    start_date=dt.date(2009,1,1)
    end_date=dt.date(2011,12,31)
    sdir='/home/pkuehn/data/TM5Inversion/xco2_mean_flexpart_region/'

    for Dtype in Dtypes:            
        path=f'/mnt/data/users/eschoema/TK5_4DVAR/Molefractions/{Dtype}/output/2009010100-2019070100/mix/'
        spath=f'{sdir}{Dtype}/'
        print(spath)
        # os.mkdir(spath)
        for date in pd.date_range(start_date, end_date):
            filepath=path+f"{date.strftime('%Y')}/{date.strftime('%m')}/mix_{date.strftime('%Y%m%d')}.nc4"
            CutAndGetPressuresTM5mix(filepath, spath, date.strftime('%Y%m%d'))
    

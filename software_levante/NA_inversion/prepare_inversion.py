import xarray as xr
import os
import pandas as pd
import datetime as dt
import numpy as np
from dateutil.relativedelta import relativedelta
from utils import get_start_date_of_week, read_grid_time_file, add_to_gosat_sounding_data, p

# save last positions as dataframe
def createLastPositionsDF(path):
    '''Get last position of each released particle, save as dataframe
    Args:
        path: path to Flexpart output directory
    Returns: nothing, saves DF_last_positions.pkl'''
    if not path[-1]=='/':
        path=path+'/'
    particle_files=[f for f in os.listdir(path) if f.startswith('partoutput_')]
    firstFile=True
    # combine data into one dataset
    print('Reading data')
    for f in particle_files:
        if firstFile:
            data=xr.open_dataset(path+f)
            firstFile=False
        else:
            temp=xr.open_dataset(path+f)
            data=xr.concat([data,temp], dim='time')
    # sort by time
    data=data.sortby('time')
    
    # get first times
    first_times = data.lon.notnull().idxmax('time')
    # select data for those fist positions, create dataframe
    data_df=data.where(data.time==first_times).to_dataframe().dropna()
    # save df
    print(f'Saving DataFrame to: {path}DF_last_positions.pkl')
    data_df.reset_index().to_pickle(path+'DF_last_positions.pkl')
    return

# calculate enhancement based on TM5-4DVar posterior fluxes transported with FLEXPART footprints
def calc_TM5flux_enhancement(footprint_dir, flux_path,gosat_path,gosat_spath='', col_name='TM5_xco2_enhancement'):
    """Calculate co2 enhancement from FLEXPART footprints
    Args:
        footprint_dir: path to directory containing FLEXPART otput
        flux_path: path to TM5-4DVar flux data containing variable total flux and coordinates
        gosat_path: path to file containing gosat sounding position, time and xco2
        gosat_spath: path to where gosat file, now including mean background, should be saved, defaults to gosat_path
        col_name: column name for enhancement data that is added to gosat file, defaults to 'TM5_xco2_enhancement
    Returns:
        xco2_enhancement_estimate: array of enhancement estimates, prior*footprint
        saves enhancement estimates into gosat sounding position file 
    """
    print(col_name)
    # read footprint data
    data=read_grid_time_file(footprint_dir)
    # get months from data
    release_months=np.unique(data.time.values.astype('datetime64[M]'))
    release_month_tuple = [(str(dt).split('-')[0], str(dt).split('-')[1]) for dt in release_months]
    release_month_tuple = [(int(year), int(month)) for year, month in release_month_tuple]
    print(f'release month: {release_month_tuple}')

    # select height
    data=data.sel(height=30).squeeze(dim=['numspec', 'nageclass'])

    # sum over time seperately for each month
    for i in range(0,len(release_month_tuple)):
        if i==0:
            data_sel=data.sel(time=f'{release_month_tuple[i][0]}-{release_month_tuple[i][1]:02d}').sum(dim='time').assign_coords(time=dt.datetime(release_month_tuple[i][0], release_month_tuple[i][1],1))
        else:
            data_temp=data.sel(time=f'{release_month_tuple[i][0]}-{release_month_tuple[i][1]:02d}').sum(dim='time').assign_coords(time=dt.datetime(release_month_tuple[i][0], release_month_tuple[i][1],1))
            data_sel=xr.concat([data_sel, data_temp], dim='time')
    # read flux
    print('reading flux data')
    flux_data=xr.open_dataset(flux_path)
    # select correct month
    flux_sel=flux_data.sel(time=[f'{year}-{month}' for (year,month) in release_month_tuple])
    enhancement=data_sel.spec001_mr*flux_sel.total_flux/data_sel.height.item()*28.9/44*1e6 # mass mixing ratio to mole fraction in ppm, *molarmass_dry_air/molarmass_co2
    enhancement=enhancement.sum(dim=['time','latitude','longitude']).values

    # add pointspec and enhancement estimate to gosat files
    add_to_gosat_sounding_data([enhancement], [col_name], gosat_path,gosat_spath)
    return enhancement
# calculate background for GOSAT meas(& get TM5_4DVar xco2 value)
def calc_TM5_background_gosat(release_dir, TM5_dir,num_parts,gosat_path,gosat_spath='',interp_method='linear', col_name='', xco2_col_name='',get_TM5_xco2_val=True, plim=100):
    """Calculate co2 background from dataframe containing last positions of all particles, path to TM5 data with pressure at boundaries 
    Args:
        last_positions_path: path to DataFrame containing, lat,lon, time and pressure of last position of all released particles of one flexpart run
        TM5_dir: path to directory where TM5 concentration data is saved, needs to include pressure at boundaries
        num_parts (int): number of particles per release, assuming one release layer per sounding position
        gosat_path: path to file containing gosat sounding position, time and xco2
        gosat_spath: path to where gosat file, now including mean background, should be saved, defaults to gosat_path
        interp_method: str for interpolation method, default: linear
        col_name: name for TM5 background column, defaults to 'TM5_background_interp_method'
        get_TM5_xco2_val: set to True if TM5-4DVar xco2 value should be saved into .csv file, default True
        xco2_col_name: name for TM5 xco2 column, defaults to 'TM5_xco2'
        plim:   upper boundary of particle releases in hPa, defults to 100hPa
    Returns:
        background_mean: array of average background for each sounding position
        saves background into gosat sounding position file
    """ 
    if col_name=='':
        col_name=f'TM5_background_{interp_method}'
    if xco2_col_name=='':
        xco2_col_name='TM5-4DVar_xco2_nearest'
    # read last positions data
    last_positions_df=pd.read_pickle(f"{release_dir}DF_last_positions.pkl")
    # check if pressure exists in df
    if not 'prs' in last_positions_df.columns:
        # caluculate height from pressure
        last_positions_df['prs']=p(last_positions_df.z)*100 # p[hPa]->p[Pa]

    # for last positions file, get timerange, read respective TM5 files
    date_min=np.min(last_positions_df.time).date()
    date_max=np.max(last_positions_df.time).date()

    # read TM5 data for that time range
    firstFile=True
    for date in pd.date_range(date_min,date_max):
        date_str=date.strftime('%Y%m%d')
        if firstFile:
            TM5_data=xr.open_dataset(f'{TM5_dir}xco2_mean_{date_str}.nc')
            firstFile=False
        else:
            temp=xr.open_dataset(f'{TM5_dir}xco2_mean_{date_str}.nc')
            TM5_data=xr.concat([TM5_data,temp], dim='times')

    # add boundaries as coordinate 
    TM5_data['boundaries']=TM5_data.boundaries
    TM5_data=TM5_data.squeeze()
    
    # xarray with only particle number as dimension
    data=last_positions_df.set_index('particle').to_xarray()
    # get value of CT data based on particle time, lat and lon using interp_method
    TM5_interp=TM5_data.interp(times=data.time, latitude=data.lat, longitude=data.lon,method=interp_method)
    # get level the particles are in based on presssure
    # idxmax gives index of first time condition is true, -1 to get level nuber from upper boundary of level
    p_levels=((TM5_interp.p_boundary<data.prs).idxmax(dim='boundaries')-1).reset_coords(names=['times','latitude','longitude'])      
    # set values of -1 to 0 (-1 if pressure lower than lowest boundary)
    p_levels['boundaries']=p_levels.boundaries.where(p_levels.boundaries != -1, 0)        # returns value from p_levels[boundary where condition is True, otherwise fills in other value (here 1)]
    # get co2 of respective level
    background=TM5_interp.mix.sel(levels=p_levels.boundaries)
    # separate different sounding positions
    background_part=background.coarsen(particle=num_parts).mean()

    # get partial xco2 above plim from TM5-4DVar data
    plim=plim *100 # hPa-> Pa
    # read data for release day
    release_day_str=release_dir[-9:-1]
    TM5_data=xr.open_dataset(f'{TM5_dir}xco2_mean_{release_day_str}.nc')
    TM5_data=TM5_data.assign_coords(boundaries=TM5_data.boundaries, levels=TM5_data.levels).squeeze()
    # determine where plim is compared to boundaries, in which level
    upper_b=((TM5_data.p_boundary<plim).idxmax(dim='boundaries')) # upper boundary of level =level+1

    # for each level weigh mixing ratio with pressure difference over that level
    # for level containing plim, get difference from plim to upper boundary
    # returns value from TM5_data.p_diff where condition is True, so all levels that are not upper_b-1, otherwise fills in plim
    TM5_data['p_diff']=TM5_data.p_diff.where(TM5_data.levels!=(upper_b-1), (plim-TM5_data.p_boundary.sel(boundaries=upper_b)))#.sel(levels=17)
    # drop lower levels
    TM5_data=TM5_data.where(TM5_data.levels>=(upper_b-1), drop=True)
    # weigh mix with p_diff, sum over levels
    TM5_data['xco2_partial']=(TM5_data.mix*TM5_data.p_diff).sum(dim='levels')/plim
    
    # get xco2 value for nearest time & position for each gosat sounding position
    # read gosat data
    gosat_data=pd.read_csv(gosat_path)
    val=[]
    TM5_xco2_val=[]
    for i in range(0,len(gosat_data)):
        temp=TM5_data.sel(latitude=gosat_data.latitude[i],longitude=gosat_data.longitude[i],times=gosat_data.time[i], method='nearest')
        psurf=temp.pressure.values[0]
        #particle contribution
        temp_p=background_part[i].item()
        # weighted sum with BG contribution from particles
        sum=temp.xco2_partial.item()*plim/psurf+ temp_p*(1-plim/psurf)
        val.append(sum)
        if get_TM5_xco2_val:
            # get xco2 value from TM5_4Dvar dataset
            TM5_xco2_val.append(temp.xco2.values[0])
    if get_TM5_xco2_val:
        add_to_gosat_sounding_data([val, TM5_xco2_val], [col_name, xco2_col_name],gosat_path,gosat_spath)
    else:
        add_to_gosat_sounding_data([val], [col_name],gosat_path,gosat_spath)  
    return val

def get_gosta_bg_xco2(start_date, end_date,flex_data_dir, gosat_dir,gosat_sdir,TM5_ds, TM5flux_path, TM5_dir, num_parts=40000):
    ''' calculate gosat background, TM5-4DVar xco2 and flux enhancement estimate using 1x1 footprint and TM5-4DVar posterior flux
    Args:
        start_date, end_date: timeperiod of measurements
        flex_data_dir: directory to flexpart output with subdirectories /YYYY_MM/Release_YYYYMMDD/
        gosat_dir: directory to gosat .csv data initially read from, with subdirectories /YYYY_MM/
        gosat_sdir: directory where gosat .csv should be saved
        TM5_ds:     TM5-4DVar dataset that should be used for background
        TM5flux_path:   path to TM5-4DVar fluxes, used for enhancement estimate with 1x1 flux*footprint
        TM5_dir:        directory containing to TM5-4DVar molefractions
        num_parts: number of particles released per measurement, defaults to 40000
    Returns:
        nothing, saves data into .csv files in specified directory
    '''
    for d in pd.date_range(start_date,end_date):
        print(f'Current date: {d.strftime("%Y%m%d")}')
        release_dir=f'{flex_data_dir}/{d.strftime("%Y_%m")}/Release_{d.strftime("%Y%m%d")}/'
        gosat_path=f"{gosat_dir}/{d.strftime('%Y_%m')}/RemoTeCv2.4.0_{d.strftime('%Y%m%d')}.csv"  
        gosat_spath=f"{gosat_sdir}/{d.strftime('%Y_%m')}/xco2_bg_{d.strftime('%Y%m%d')}.csv"                         
        if not os.path.isdir(f"{gosat_sdir}/{d.strftime('%Y_%m')}"):
            os.makedirs(f"{gosat_sdir}/{d.strftime('%Y_%m')}")
        # check if there are GOSAT measurements for that day
        if os.path.isfile(gosat_path):
            # get last positions of particles
            last_positions_path=f"{release_dir}/DF_last_positions.pkl"
            # check if DF_last_positions exists, create if doesnt exist
            if not os.path.isfile(last_positions_path):
                print('create DF_last_positions.pkl')
                createLastPositionsDF(release_dir)
            
            # calculate co2 enhancement
            # check if .csv already exists at desired saving location
            if os.path.isfile(gosat_spath): # read csv from gosat_spath
                calc_TM5flux_enhancement(release_dir, TM5flux_path,gosat_spath,gosat_spath, col_name=f'TM5_{TM5_ds}_enhancement')
            else:    # read csv from gosat_path
                calc_TM5flux_enhancement(release_dir, TM5flux_path,gosat_path,gosat_spath, col_name=f'TM5_{TM5_ds}_enhancement')
            
            
            # calc 100 hPa background & get TM5-4DVar xco2 
            calc_TM5_background_gosat(release_dir, TM5_dir,num_parts,gosat_spath, gosat_spath, col_name=f'TM5_{TM5_ds}_background',xco2_col_name=f'TM5_{TM5_ds}_xco2')

def calc_interpolated_TM5_meas_values(start_date, end_date, TM5_dir, gosat_dir, is_dir, bg_str='RemoTeC_2.4.0+IS'):
    ''' Calculate interpolated molefraction values from TM5-4DVar
    Args:
        start_date (datetime.date): start date for releases that should be used
        end_date (datetime.date): end date for releases that should be used (including this day)
        TM5_dir: path to TM5-4DVar molefractions data
        gosat_dir: path to Gosat csv directories 
        is_dir: path to insitu csv directories 
        bg_str: string indicating which TM5-4DVar dataset should be used, defaults to 'RemoTeC_2.4.0+IS'
    Returns: nothing, saves interpolated TM5-4DVar data into csv files
    '''
    # Function to interpolate a single row
    def interpolate_xco2_point(row):
        val = TM5_data.xco2.interp(
            times=row['time'],
            latitude=row['latitude'],
            longitude=row['longitude'],
            method='linear'  # or 'nearest' if needed
        )
        return val.values.item()  # extract scalar
    def interpolate_co2_point(row):
        temp = TM5_data.mix.interp(
            times=row['time'],
            latitude=row['latitude'],
            longitude=row['longitude'],
            method='linear'  # or 'nearest' if needed
        )
        val = temp.swap_dims({'levels':'p_level'}).interp(p_level=row['p_intake_height[hPa]'])
        return val.values.item()  # extract scalar
    
    for date in pd.date_range(start_date, end_date):
        # read data
        TM5_data=xr.open_dataset(f'{TM5_dir}xco2_mean_{date.strftime("%Y%m%d")}.nc').squeeze()
        # gosat
        gosat_path=f'{gosat_dir}/{date.strftime("%Y_%m")}/xco2_bg_{date.strftime("%Y%m%d")}.csv'
        if os.path.isfile(gosat_path):  # check if file exists, there are some days without measurements
            gosat_csv=pd.read_csv(gosat_path, parse_dates=['time'])
            # Apply interpolation function to each row in the DataFrame
            gosat_csv[f'TM5_{bg_str}_xco2_interpolated'] = gosat_csv.apply(interpolate_xco2_point, axis=1)
            # save csv files
            gosat_csv.to_csv(gosat_path, index=None)
        # insitu
        else: 
            print(f'no gosat file for {date.strftime("%Y%m%d")}')
        is_path=f'{is_dir}/{date.strftime("%Y_%m")}/co2_bg_{date.strftime("%Y%m%d")}.csv'
        if os.path.isfile(is_path):
            is_csv=pd.read_csv(is_path, parse_dates=['time'])
            # get mid pressure for level from pressure at boundaries
            # /100 to get pressure in hPa
            TM5_data['p_level']=((TM5_data.p_boundary.sel(boundaries=slice(0,25))+TM5_data.p_boundary.sel(boundaries=slice(1,26)))/2/100).swap_dims({'boundaries':'levels'})
            TM5_data=TM5_data.assign_coords(p_level=TM5_data.p_level)
            is_csv[f'TM5_{bg_str}_co2_interpolated'] = is_csv.apply(interpolate_co2_point, axis=1)
            is_csv.to_csv(is_path, index=None)
            print(f'saved for {date.strftime("%Y%m%d")}')
        else:
            print(f'no is file for {date.strftime("%Y%m%d")}')
    return

def calc_TM5_background_for_ISmeas(release_dir, TM5_dir,num_parts,is_path,is_spath='',interp_method='linear', col_name='', co2_col_name=''):
    """Calculate co2 background for in-situ measurements from dataframe containing last positions of all particles, path to TM5 data with pressure at boundaries 
    Args:
        last_positions_path: path to DataFrame containing, lat,lon, time and pressure of last position of all released particles of one flexpart run
        TM5_dir: path to directory where TM5 concentration data is saved, needs to include pressure at boundaries
        num_parts (int): number of particles per release, assuming one release layer per sounding position
        is_path: path to file containing insitu measurement position, time and co2
        is_spath: path to where insitu file, now including mean background, should be saved, defaults to is_path
        interp_method: str for interpolation method, default: linear
        col_name: name for TM5 background column, defaults to 'TM5_background_interp_method'
        get_TM5_xco2_val: set to True if TM5-4DVar xco2 value should be saved into .csv file, default True
        xco2_col_name: name for TM5 xco2 column, defaults to 'TM5_xco2'
    Returns: nothing
        saves background into insitu measurement position file
    """ 
    if col_name=='':
        col_name=f'TM5_background_{interp_method}'
    if co2_col_name=='':
        co2_col_name='TM5-4DVar_co2_nearest'
    # read last positions data
    last_positions_df=pd.read_pickle(f"{release_dir}DF_last_positions.pkl")
    # check if pressure exists in df
    if not 'prs' in last_positions_df.columns:
        # caluculate height from pressure
        last_positions_df['prs']=p(last_positions_df.z)*100 # p[hPa]->p[Pa]
    
    # for last positions file, get timerange, read respective TM5 files
    date_min=np.min(last_positions_df.time).date()
    date_max=np.max(last_positions_df.time).date()

    # read TM5 data for that time range
    TM5_list=[]
    for date in pd.date_range(date_min,date_max):
        date_str=date.strftime('%Y%m%d')
        TM5_data=xr.open_dataset(f'{TM5_dir}xco2_mean_{date_str}.nc')
        TM5_list.append(TM5_data)
    TM5_data=xr.concat(TM5_list, dim='times')
    # add boundaries as coordinate 
    TM5_data['boundaries']=TM5_data.boundaries
    TM5_data=TM5_data.squeeze()

    # xarray with only particle number as dimension
    data=last_positions_df.set_index('particle').to_xarray()
    # get value of CT data based on particle time, lat and lon using interp_method
    TM5_interp=TM5_data.interp(times=data.time, latitude=data.lat, longitude=data.lon,method=interp_method)
    # get level the particles are in based on presssure
    # idxmax gives index of first time condition is true, -1 to get level nuber from upper boundary of level
    p_levels=((TM5_interp.p_boundary<data.prs).idxmax(dim='boundaries')-1).reset_coords(names=['times','latitude','longitude'])      
    # set values of -1 to 0 (-1 if pressure lower than lowest boundary)
    p_levels['boundaries']=p_levels.boundaries.where(p_levels.boundaries != -1, 0)        # returns value from p_levels[boundary where condition is True, otherwise fills in other value (here 1)]
    # get co2 of respective level
    background=TM5_interp.mix.sel(levels=p_levels.boundaries)
    # separate different sounding positions
    background_part=background.coarsen(particle=num_parts).mean()

    # get co2 value for nearest time & position for each insitu measurement position
    # read insitu data
    is_data=pd.read_csv(is_path)        # , index_col=0
    # add background values to dataframe
    is_data[col_name]=background_part.values

    TM5_co2_val=[]
    p_intake_height=[]
    for i in range(0,len(is_data)):
        temp=TM5_data.sel(latitude=is_data.latitude[i],longitude=is_data.longitude[i],times=is_data.time[i], method='nearest')
        # surface pressure
        psurf=temp.pressure.item()
        # get pressure at intake_height
        p_intake=p(is_data['intake_height[magl]'][i], p0=psurf)
        p_intake_height.append(p_intake)
        # get co2 val from corresponting layer from temp_TM5_data
        # determine level the of the intake_height
        p_level=((temp.p_boundary<p_intake).idxmax(dim='boundaries')-1).item() # idxmax gives index of first time condition is true, -1 to get level number from upper boundary of level
        # print(p_level)
        # get corresponding mixing ratio
        TM5_co2_val.append(temp.mix.sel(levels=p_level).item())
    # add to is_data df
    is_data['p_intake_height[hPa]']=np.array(p_intake_height)/100
    is_data[co2_col_name]=TM5_co2_val
    print(f'saving to {is_spath}')
    is_data.to_csv(is_spath, index=None)
    return

def process_flexpart_runs(start_date, end_date,flex_dir,gosat_dir,is_dir, TM5flux_dir,TM5_molefrac_dir, num_parts=40000,plim=100, ds_list=['RemoTeC_2.4.0+IS']):
    ''' calculate background, TM5-4DVar values for gosat and insitu flexpart runs, get last positions and calculate weekly footprint sums
    Args:
        start_date, end_date: dt.date() defining the time period
        flex_dir: path to Flexpart output directories, with subdirectories 'insitu/yyyy_mm/Release_yyyymmdd', 'RemoTeCv240/...'
        gosat_dir: path to gosat measurement .csv files, with subdirectories yyyy_mm
        is_dir: path to insitu measurement .csv files, with subdirectories yyyy_mm
        TM5flux_dir: path to TM5-4DVar flux directory
        TM5_molefrac_dir: path to TM5-4DVar molefractions
        
        optional:
        num_parts: number of particles per release, defaults so 40000
        plim:   upper boundary of particle releases in hPa, defults to 100hPa
        ds_list: list of TM5-4DVar reference fluxes that should be used for background calculations, defaults to 'RemoTeC_2.4.0+IS', other ones are not relevant at the moment
                    chose dataset from ['RemoTeC_2.4.0+IS', 'ACOS+IS','IS', 'RemoTeC_2.4.0']
    Returns:
        Nothing, will create insitu/TM5-4DVar_estimate/yyyy_mm and RemoTeCv240/TM5-4DVar_estimate/yyyy_mm subdirectories to save .csv files
    '''
    # get last_positions dataframe, calculate background, enhancement from TM5-4DVar posterior 1x1 fluxes * 1x1 footprint
    for ds in ds_list:
        TM5flux_path=f'{TM5flux_dir}/flux_1x1_{ds}_cut.nc'
        TM5_dir=f"{TM5_molefrac_dir}/{ds}/"
        
        # gosat measurements:
        # calculate background, TM5-4DVAR xco2 and enhancement from footprints an 1x1 fluxes
        flex_data_dir=f'{flex_dir}/RemoTeCv240/'
        gosat_sdir=f"{flex_dir}/RemoTeCv240/TM5-4DVar_estimate/"
        get_gosta_bg_xco2(start_date, end_date,flex_data_dir, gosat_dir,gosat_sdir,ds, TM5flux_path, TM5_dir)
    
        # insitu measurements:
        # calculate background
        is_sdir=f'{flex_dir}/insitu/TM5-4DVar_estimate'
        for date in pd.date_range(start_date, end_date):
            release_dir=f'{flex_dir}/insitu/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/'
            is_path=f'{is_dir}/{date.strftime("%Y_%m")}/ISpositions_{date.strftime("%Y%m%d")}.csv'
            if not os.path.isdir(f'{is_sdir}/{date.strftime("%Y_%m")}'):
                os.makedirs(f'{is_sdir}/{date.strftime("%Y_%m")}')
            is_spath=f'{is_sdir}/{date.strftime("%Y_%m")}/co2_bg_{date.strftime("%Y%m%d")}.csv'
            
            last_positions_path=f"{release_dir}/DF_last_positions.pkl"
            # check if DF_last_positions exists, create if doesnt exist
            if not os.path.isfile(last_positions_path):
                print('create DF_last_positions.pkl')
                createLastPositionsDF(release_dir)
                # get background for insitu measurements
            calc_TM5_background_for_ISmeas(release_dir, TM5_dir,num_parts,is_path,is_spath, col_name=f'TM5_{ds}_background',co2_col_name=f'TM5_{ds}_co2')
        # get interpolated TM5-4DVar molefractions
        calc_interpolated_TM5_meas_values(start_date, end_date, TM5_dir, gosat_sdir, is_sdir, bg_str)

def get_frac_remaining_particles(dir_path, csv_dir,csv_sdir,file_str, start_date, end_date, num_parts=40000):
    ''' Calculate fraction of remaining particles from Flexpart runs for each release, 
            saves them to .csv files with measurement data and creates .nc file with fraction for all releases
    Args:
        dir_path: path to flexpart release directories with subdirectories YYYY_MM/Release_YYYYMMDD
        csv_dir: path to csv files with subdirectories YYYY_MM
        csv_sdir: path to where csv files should be saved, creates YYYY_MM subdirectories if they dont exist
        file_str: file name, eg. xco2_bg for gosat files
        start_date, end_date: dt.date objects defining the time period
        num_parts: number of particles per release, defaults to 40000
    Returns: nothing, saves csv files and all fractions to one .nc file
    '''
    data_list=[]
    for date in pd.date_range(start_date, end_date):
        # check that release dir exists
        release_path=f'{dir_path}/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/'
        if os.path.isdir(release_path):
            # get min partoutput file
            partoutput_files=[f for f in os.listdir(release_path) if f.startswith('partoutput') if f.endswith('.nc')]
            partoutput_files.sort()
            # read data from min partoutput
            partout_data=xr.open_dataset(release_path+partoutput_files[0])
            partout_data=partout_data.sel(time=partout_data.time.min(), drop=True)
            # unstack coordinates to particles and release number
            release_num=(partout_data.particle.values - 1) // num_parts
            part=(partout_data.particle.values - 1) % num_parts + 1
            partout_data=partout_data.assign_coords(release_num=("particle", release_num))
            partout_data = partout_data.assign_coords(part=("particle", part))
            partout_data=partout_data.swap_dims({"particle": "part"}).drop_vars('particle')
            partout_data=partout_data.set_index(particle_temp=["part", "release_num"]).unstack("particle_temp")
            # count remaining particles
            partout_data['particles_remaining']=(partout_data.z.count(dim='part')).assign_attrs(description='number of particles remaining at end of simulation')
            # in percent
            partout_data['remaining']=(partout_data.particles_remaining/num_parts).assign_attrs(description='% of particles remaining at end of simulation')
            # add that to .csv file
            file_path=f'{csv_dir}/{date.strftime("%Y_%m")}/{file_str}_{date.strftime("%Y%m%d")}.csv'
            data=pd.read_csv(file_path, index_col=0)
            data['frac_remaining']=partout_data['remaining'].values
            if not os.path.isdir(f'{csv_sdir}/{date.strftime("%Y_%m")}'):
                os.makedirs(f'{csv_sdir}/{date.strftime("%Y_%m")}')
            data.to_csv(f'{csv_sdir}/{date.strftime("%Y_%m")}/{file_str}_{date.strftime("%Y%m%d")}.csv')

            # without time as dimension
            data_list.append(partout_data.remaining)
    if len(data_list)==0:
        print('no data found, problem with given directory?')
        print(release_path)
    else:
        data=xr.concat(data_list, dim='time')
        # save to nc file
        spath=f'{dir_path}/remaining_particles_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.nc'
        print(f'saving to {spath}')
        data.to_netcdf(spath)

def cut_total_tm5_4DVar_flux(ds_path, ):
    ''' Cut TM5-4DVar data to desired region, calculate total flux, save dataset
    Args:
        ds_path: path to the dataset
        region= list of [lat_min, lat_max, lon_min, lon_max]
    Returns: nothing, saves cut dataframe
    '''
    flux_data=xr.open_dataset(ds_path)


    # calculate total flux
    flux_data['total_flux']=flux_data.total_flux.assign_attrs(units='kgCO2/(m^2 s)')
    flux_data.to_netcdf(f'{ds_path[:-3]}_cut.nc')
    return

def prep_TM5_4DVar_flux(flux_path, region, res):
    ''' Cut TM5-4DVar data to desired region, get total flux, adjusts flux units as needed for the inversion and calculate average flux for a coarser spatial grid
    Args:
        flux_path: path to where flux data is saved
        res: desired resolution, e.g. 2 for 2x2 spatial grid
    Returns:
        saves netcdf file to flux_path with 'flux_1x1'  replaced by ne spatial grid resoltion 'flux_2x2'
    '''
    data=xr.open_dataset(flux_path)
    # shift coordinates from 0-360 and 0-180 to -180 - 180 and -90 - 90
    data = data.assign_coords(time=("months",[dt.datetime(year, month,1) for year, month in data.month_tuple.values]), latitude=data.latitude-89.5, longitude=data.longitude-179.5).swap_dims({"months": "time"})
    # cut flexpart region
    data=data.sel(latitude=slice(region[0],region[1]),longitude=slice(region[2],region[3]))
    # get total flux
    data['total_flux']=(data.CO2_flux_nee+data.CO2_flux_fire+data.CO2_flux_oce+data.CO2_flux_fos)
    
    data=data[['grid_cell_area','CO2_flux_nee','CO2_flux_fire','CO2_flux_oce','CO2_flux_fos','total_flux']]
    # adjust units for flux components
    cols=['CO2_flux_nee','CO2_flux_fire','CO2_flux_oce','CO2_flux_fos','total_flux']
    for col in cols:
        data[col]=(data[col]*44/12*1e-3/(24*60*60)).assign_attrs(units='kgCO2/(m^2 s)')

    # area weighted mean over coarser grid
    data[cols]=data[cols]*data.grid_cell_area
    data=data.coarsen(latitude=res, longitude=res).sum().squeeze()
    data[cols]=(data[cols]/data.grid_cell_area)
    for col in cols:
        data[col]=data[col].assign_attrs(units='kgCO2/(m^2 s)')
    # save dataset
    spath=f"{flux_path.replace('_1x1_', f'_{res}x{res}_')[:-3]}_cut.nc"
    print(f'save preprocessed TM5-4DVar fluxes to {spath}')
    data.to_netcdf(spath)
    return

# get weekly priors
def get_weekly_TM5_4DVarflux(ds_path):
    ''' get weekly mean from TM5-4DVar flux dataset
    Args: 
        ds_path: path to flux dataset
    Returns:
        nothing, saves weekly average as dataset
    '''
    # read TM5-4DVar flux
    ds=xr.open_dataset(ds_path)
    # timerange of the TM5-4DVar dataset
    TM5_start_date=dt.date(2009,1,1)
    TM5_end_date=dt.date(2019,6,30)
    # Convert monthly data to daily by forward-filling values
    ds = ds.reindex(time=pd.date_range(TM5_start_date,TM5_end_date), method='ffill')

    # Create the 7-day period bins
    period_bins = pd.date_range(start=get_start_date_of_week(TM5_start_date), end=get_start_date_of_week(TM5_end_date)+dt.timedelta(days=7), freq="7D")
    # Cut the time array into the defined 7-day bins
    time_bins = pd.cut(ds.time, bins=period_bins, right=False, labels=period_bins[:-1],include_lowest=True)
    # Assign the new time bins as coordinates
    ds=ds.assign_coords(time=time_bins)

    # weekly mean, only different value for weeks across two months
    # the timestamp corresponds to the first day of the week
    ds_mean=ds.groupby("time").mean('time')
    ds_mean['time']=ds_mean.time.assign_attrs(description='start date of the week that is averaged over')
    spath=ds_path.replace('.nc', '_weekly.nc')
    print(f'saving data to {spath}')
    ds_mean.to_netcdf(spath)
    return

# main
if __name__ == "__main__":  
    start_date, end_date =dt.date(2010,6,1), dt.date(2010,6,6)    # time period of measurements
    # data paths
    # TM5-4DVar molefractions with selected dataset as subdirectory
    TM5_molefrac_dir='/work/bb1170/RUN/b382762/data/TM5Inversion/co2_xco2_mean/'
    bg_str = 'RemoTeC_2.4.0+IS'         # TM5-4DVar dataset, chose from ['RemoTeC_2.4.0+IS', 'ACOS+IS','IS', 'prior]
    TM5_dir=f'{TM5_molefrac_dir}{bg_str}/'
    # directory with TM5-4DVar fluxes
    TM5flux_dir='/work/bb1170/RUN/b382762/data/FarewellPackage_test/TM5-4DVar/' 
    # path to high resolution (interpolated) scaling data
    scaling_data_path='/work/bb1170/RUN/b382762/data/TM5Inversion/high_res_total_scaling_RemoTeC+IS.nc'

    # FLEXPART
    flex_dir='/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/' # path to directory with RemoTeC and insitu Flexpart release subdirectories
    num_parts=40000     # number of particles per measurement

    # measurement data paths
    # path to soundingposition.csv files
    gosat_dir='/work/bb1170/RUN/b382762/data/FarewellPackage_test/GOSAT/'      
    is_dir=f'/work/bb1170/RUN/b382762/data/FarewellPackage_test/ObsPack/'
    # path to where .csv files should be saved, will creade subdirectories YYYY_mm/co2_bg_YYYYMMDD.csv
    # should not be changed
    gosat_csv_sdir=f'{flex_dir}/RemoTeCv240/TM5-4DVar_estimate/'
    is_csv_dir=f'{flex_dir}/insitu/TM5-4DVar_estimate/'   
    # if change folder substructure/ folder names, need to adapt that in all functions and sections below
    # e.g flex_dir/insitu/YYYY_MM/..., RemoTeCv240, or TM5-4DVar_estimate

    # select funtions to run
    PREP_TM5_4DVAR_REF_FLUXES=True
    
    GET_IS_BG_CO2=False # get background and TM5-4DVar moelfractions, save as .csv files in specified directory
    GET_XCO2_VALS=False
    GET_INTERP_TM5_VALS=False
    PROCESS_FLEX_RUNS=True  # combines the three steps above
    GET_FRAC_REM=True
    
    GET_HIGH_RES_FOOTPRINTS=True
    GET_TM5_4DVAR_SCALED_FOOTPRINTS=True    # will also calculate unscaled weekly footprints, optionally: remove hourly footprint files for storage reasons
    GET_WEEKLY_NO_SCALING=True         # use if no scaling should be applied to the footprints, adapt path in COARSEN_HIGH_RES_FOOTPRINT in that case
    COARSEN_HIGH_RES_FOOTPRINT=True     # if True, adapt dir_path
    
    if PREP_TM5_4DVAR_REF_FLUXES:
        # cuts specified region, coarsens to desired resolution, saves weekly files
        # for dataset selected for background calculation and prior
        for temp in [bg_str, 'prior']:
            flux_path=f'{TM5flux_dir}/flux_1x1_{temp}.nc'
            region=[6,66,-146,-56] # somewhat larger than Flexpart region 12,56,-134,-62
            res=2
            # cut fluxes to desired region for 1x1 and res resolution
            prep_TM5_4DVar_flux(flux_path, region, 1)
            prep_TM5_4DVar_flux(flux_path, region, res)
            # get weekly fluxes for 1x1 and 2x2
            get_weekly_TM5_4DVarflux(f"{flux_path[:-3]}_cut.nc")
            get_weekly_TM5_4DVarflux(f"{flux_path.replace('_1x1_', f'_{res}x{res}_')[:-3]}_cut.nc")
    
    # calculate DF_last_positions.pkl, get background, xco2 enhancements by selecting nearest and by interpolating
    if GET_IS_BG_CO2:
        for date in pd.date_range(start_date, end_date):
            release_dir=f'{flex_dir}/insitu/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/'
            is_path=f'{is_dir}{date.strftime("%Y_%m")}/ISpositions_{date.strftime("%Y%m%d")}.csv'
            is_spath=f'{is_csv_dir}{date.strftime("%Y_%m")}/co2_bg_{date.strftime("%Y%m%d")}.csv'
            if not os.path.isdir(f'{is_csv_dir}{date.strftime("%Y_%m")}'): # create monthly subdirectories if they dont exist
                os.makedirs(f'{is_csv_dir}{date.strftime("%Y_%m")}')
            if not os.path.isfile(f'{release_dir}DF_last_positions.pkl'): # create last_positions file if it doesnt exists
                print('create DF_last_positions.pkl')
                createLastPositionsDF(release_dir)
            calc_TM5_background_for_ISmeas(release_dir, TM5_dir,num_parts,is_path,is_spath, col_name=f'TM5_{bg_str}_background',co2_col_name=f'TM5_{bg_str}_co2')
    if GET_XCO2_VALS:
        # path to TM5-4DVar fluxes, used for enhancement estimate with 1x1 flux*footprint
        TM5flux_path_1x1=f'{TM5flux_dir}/flux_1x1_{bg_str}_cut.nc'    
        get_gosta_bg_xco2(start_date, end_date,f'{flex_dir}/RemoTeCv240', gosat_dir,gosat_csv_sdir,bg_str, TM5flux_path_1x1, TM5_dir)
    if GET_INTERP_TM5_VALS:
        calc_interpolated_TM5_meas_values(start_date, end_date, TM5_dir, gosat_csv_sdir, is_csv_dir, bg_str)
    
    # three steps above combined in this one
    if PROCESS_FLEX_RUNS:
        process_flexpart_runs(start_date, end_date,flex_dir,gosat_dir,is_dir, TM5flux_dir,TM5_molefrac_dir)

    if GET_FRAC_REM:
        # insitu releases
        file_str='co2_bg'   # .csv file name for insitu measurements
        get_frac_remaining_particles(f'{flex_dir}/insitu/', is_csv_dir,is_csv_dir,file_str, start_date, end_date)

        # gosat releases
        file_str='xco2_bg' # .csv file name for gosat measurements
        get_frac_remaining_particles(f'{flex_dir}/RemoTeCv240/', gosat_csv_sdir,gosat_csv_sdir,file_str, start_date, end_date)

    # get 1x1 hourly footprints with measurement data into one ds per month, needed for diurnal cacle
    if GET_HIGH_RES_FOOTPRINTS:
        for month_start in pd.date_range(start_date, end_date, freq='MS'):
            month_end= month_start+relativedelta(months=1, days=-1)
            print(month_start, month_end)
            dir_list=[f'{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/' for date in pd.date_range(month_start, month_end) if date in pd.date_range(month_start, month_end)]
            
            # for insitu data
            flexpart_is_output_path= f'{flex_dir}/insitu/'
            is_files=[flexpart_is_output_path+d+f for d in dir_list if os.path.isdir(flexpart_is_output_path+d) for f in os.listdir(flexpart_is_output_path+d) if f.startswith('grid_time_')]
            ds=[]
            for f in is_files:
                temp=xr.open_dataset(f).sel(height=30)[['spec001_mr']].squeeze(dim='nageclass')
                temp['release_num']=(('pointspec'), temp.pointspec.values)
                temp['release_day']=pd.to_datetime(f[-36:-28])
                ds.append(temp)
            is_footprint=xr.concat(ds, dim='pointspec')
            is_footprint['pointspec']=is_footprint.pointspec.values
            # add meas info to this dataset
            data_list=[]
            for date in pd.date_range(month_start, month_end):    
                # check if file exists for that date
                path=f'{flexpart_is_output_path}/TM5-4DVar_estimate/{date.strftime("%Y_%m")}/co2_bg_{date.strftime("%Y%m%d")}.csv'
                if os.path.isfile(path):
                    data=pd.read_csv(path, parse_dates=['time'])
                    data_list.append(data)
            is_data=pd.concat(data_list)
            is_data=is_data.reset_index(names='release_num')
            # check that release dates match
            if (is_footprint.release_day.dt.date.values==is_data.time.dt.date).all():
                # check that release numbers match 
                if (is_footprint.release_num.values==is_data.release_num).all():
                    # rename columns
                    is_data.rename(columns={'latitude':'release_lat', 'longitude':'release_lon', 'time':'release_time'}, inplace=True)
                    is_data.columns
                    # add release data to footprints
                    cols=['release_time', 'file', 'release_lat', 'release_lon','co2_val[ppm]', 'elevation[masl]', 'intake_height[magl]',
                        'TM5_RemoTeC_2.4.0+IS_background', 'p_intake_height[hPa]','TM5_RemoTeC_2.4.0+IS_co2', 'frac_remaining']
                    for col in cols:
                        is_footprint[col]=(('pointspec'), is_data[col])
                    if not os.path.isdir(f'{flexpart_is_output_path}prep_footprints/hourly'):
                        os.makedirs(f'{flexpart_is_output_path}prep_footprints/hourly')
                    s_file_path=f'{flexpart_is_output_path}prep_footprints/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                    is_footprint.to_netcdf(s_file_path)
                    print(f'saved to {s_file_path}')
                else:
                    print('Problem with Release numbers')
            else: 
                print('Problem with Release days')
        
            # for gosat data
            flexpart_gosat_output_path= f'{flex_dir}/RemoTeCv240/'
            gosat_files=[flexpart_gosat_output_path+d+f for d in dir_list if os.path.isdir(flexpart_gosat_output_path+d) for f in os.listdir(flexpart_gosat_output_path+d) if f.startswith('grid_time_')]
            ds=[]
            for f in gosat_files:
                temp=xr.open_dataset(f).sel(height=30)[['spec001_mr']].squeeze(dim='nageclass')
                temp['release_num']=(('pointspec'), temp.pointspec.values)
                temp['release_day']=pd.to_datetime(f[-36:-28])
                ds.append(temp)
            gosat_footprint=xr.concat(ds, dim='pointspec')
            gosat_footprint['pointspec']=gosat_footprint.pointspec.values
            # add meas info to this dataset
            data_list=[]
            for date in pd.date_range(month_start, month_end):    
                # check if file exists for that date
                path=f'{flexpart_gosat_output_path}/TM5-4DVar_estimate/{date.strftime("%Y_%m")}/xco2_bg_{date.strftime("%Y%m%d")}.csv'
                if os.path.isfile(path):
                    data=pd.read_csv(path, parse_dates=['time'])
                    data_list.append(data)
            gosat_data=pd.concat(data_list)
            gosat_data=gosat_data.reset_index(names='release_num')
            # check that release dates match
            if (gosat_footprint.release_day.dt.date.values==gosat_data.time.dt.date).all():
                # check that release numbers match 
                if (gosat_footprint.release_num.values==gosat_data.release_num).all():
                    # rename columns
                    gosat_data.rename(columns={'latitude':'release_lat', 'longitude':'release_lon', 'time':'release_time'}, inplace=True)
                    gosat_data.columns
                    # add release data to footprints
                    cols=['release_time', 'release_lat', 'release_lon','xco2', 'xco2_err',
                        'TM5_RemoTeC_2.4.0+IS_background', 'TM5_RemoTeC_2.4.0+IS_xco2', 'frac_remaining']
                    for col in cols:
                        gosat_footprint[col]=(('pointspec'), gosat_data[col])
                    if not os.path.isdir(f'{flexpart_gosat_output_path}prep_footprints/hourly/'):
                        os.makedirs(f'{flexpart_gosat_output_path}prep_footprints/hourly/')
                    s_file_path=f'{flexpart_gosat_output_path}prep_footprints/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                    gosat_footprint.to_netcdf(s_file_path)
                    print(f'saved to {s_file_path}')
                else:
                    print('Problem with Release numbers')
            else: 
                print('Problem with Release days')
    # apply TM5-4DVar diurnal cycle scaling, add to weekly
    if GET_TM5_4DVAR_SCALED_FOOTPRINTS:
        print('reading scaling data')
        scaling_data=xr.open_dataset(scaling_data_path)
        for dir_str in ['insitu', 'RemoTeCv240']:  #
            for month_start in pd.date_range(start_date, end_date, freq='MS'):
                month_end= month_start+relativedelta(months=1, days=-1)
                # read footprint data
                data_path=f'{flex_dir}/{dir_str}/prep_footprints/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                print('reading footprint data')
                data=xr.open_dataset(data_path)
                # multiply with TM5-4DVar total scaling factors from prior
                data['spec001_mr_scaled']=(data['spec001_mr']*scaling_data.scaling_total).assign_attrs(description='Flexpart footprint in units necessary for inversion, scaled with total diurnal scaling factor from TM5-4DVar prior fluxes (sc_tot = bio_h/tot_mean + rest_mean/tot_mean)', units='s m^2/kg')
    
                # save dataset
                data_spath=f'{flex_dir}/{dir_str}/prep_footprints/scaled_hourly/high_res_scaled_footprints_{month_start.strftime("%Y_%m")}.nc'
                if not os.path.isdir(f'{flex_dir}/{dir_str}/prep_footprints/scaled_hourly'):
                    os.makedirs(f'{flex_dir}/{dir_str}/prep_footprints/scaled_hourly')
                print(f'saving preprocessed hourly footprints to {data_spath}')
                data.to_netcdf(data_spath, mode="w")
                print('saving succesfull')

                # get weekly sum
                # Create the 7-day period bins
                period_bins = pd.date_range(start=get_start_date_of_week(month_start-dt.timedelta(days=10)), end=(get_start_date_of_week(month_end)+dt.timedelta(days=7)), freq="7D")
                # Cut the time array into the defined 7-day bins
                time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=period_bins[:-1])
                # Assign the new time bins as coordinates
                data=data.assign_coords(time=time_bins)
                # weekly sum
                # the timestamp corresponds to the first day of the week
                data=data.groupby("time").sum('time')
                # save dataset
                spath=f"{flex_dir}/{dir_str}/prep_footprints/scaled_weekly/high_res_scaled_footprints_{month_start.strftime('%Y_%m')}_weekly.nc"
                if not os.path.isdir(f'{flex_dir}/{dir_str}/prep_footprints/scaled_weekly'):
                    os.makedirs(f'{flex_dir}/{dir_str}/prep_footprints/scaled_weekly')
                print(f'saving preprocessed weekly footprints to {spath}')
                data.to_netcdf(spath, mode="w")
                print('saving succesfull')
                
                del data
                
                # optional, to limit storage space
                # remove old high_res scaling file
                # print(f'deleting old hourly footprint file: {data_path}')
                # os.remove(data_path)
    if GET_WEEKLY_NO_SCALING:
        for dir_str in ['insitu', 'RemoTeCv240']:  #
            for month_start in pd.date_range(start_date, end_date, freq='MS'):
                month_end= month_start+relativedelta(months=1, days=-1)
                # read footprint data
                data_path=f'{flex_dir}/{dir_str}/prep_footprints/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                print('reading footprint data')
                data=xr.open_dataset(data_path)

                # get weekly sum
                # Create the 7-day period bins
                period_bins = pd.date_range(start=get_start_date_of_week(month_start-dt.timedelta(days=10)), end=(get_start_date_of_week(month_end)+dt.timedelta(days=7)), freq="7D")
                # Cut the time array into the defined 7-day bins
                time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=period_bins[:-1])
                # Assign the new time bins as coordinates
                data=data.assign_coords(time=time_bins)
                # weekly sum
                # the timestamp corresponds to the first day of the week
                data=data.groupby("time").sum('time')
                # save dataset
                spath=f"{flex_dir}/{dir_str}/prep_footprints/weekly/high_res_footprints_{month_start.strftime('%Y_%m')}_weekly.nc"
                if not os.path.isdir(f'{flex_dir}/{dir_str}/prep_footprints/weekly'):
                    os.makedirs(f'{flex_dir}/{dir_str}/prep_footprints/weekly')
                print(f'saving preprocessed weekly footprints to {spath}')
                data.to_netcdf(spath, mode="w")
                print('saving succesfull')
                del data
    # combine into one ds for entire time period,  coarsen to 2x2 and 4x4 resolution
    if COARSEN_HIGH_RES_FOOTPRINT:
        is_cols=['release_num','release_day','release_time','file','release_lat','release_lon','co2_val[ppm]','elevation[masl]',
                 'intake_height[magl]','TM5_RemoTeC_2.4.0+IS_background','p_intake_height[hPa]','TM5_RemoTeC_2.4.0+IS_co2','frac_remaining']
        gosat_cols=['release_num', 'release_lat', 'release_lon', 'xco2', 'xco2_err', 'TM5_RemoTeC_2.4.0+IS_background', 'TM5_RemoTeC_2.4.0+IS_xco2', 'frac_remaining']
        dir_str=['insitu', 'RemoTeCv240']
        col_list=[is_cols, gosat_cols]
        
        # for both is and gosat meas
        for i in range(0,len(dir_str)):
            dir_path=f'{flex_dir}/{dir_str[i]}/prep_footprints/scaled_weekly'   
            # dir_path=f'{flex_dir}/{dir_str[i]}/prep_footprints/weekly'   
            cols=col_list[i]
            
            # Preprocess function: only select first timestep for specific cols
            def preprocess_get_first_timestep(ds):
                # Drop unwanted variables early
                for var in cols:
                    # if var in ds:
                    ds[var] = ds[var].isel(time=0, drop=False)  # keep time dimension (size 1)
                return ds
            
            # Build the list of filepaths
            file_list = [f"{dir_path}/high_res_scaled_footprints_{month_start.strftime('%Y_%m')}_weekly.nc" 
                         for month_start in pd.date_range(start_date, end_date, freq='MS')]
            # sort file list
            file_list.sort()
            
            # combine into one ds
            # Now use open_mfdataset
            data = xr.open_mfdataset(
                file_list,
                combine="nested",
                concat_dim="pointspec",  # Stack across pointspec
                preprocess=preprocess_get_first_timestep,             
                chunks="auto",  # Open lazily with dask
                join="outer",   # Align time, lat, lon if needed
            )
            
            # Fill NaNs only for numeric variables
            numeric_vars = [v for v in data.data_vars if np.issubdtype(data[v].dtype, np.number)]
            data[numeric_vars] = data[numeric_vars].fillna(0)
            print(f'read data for {dir_str[i]}')
            
            data.to_netcdf(f"{dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_weekly.nc")
            print(f"saved {dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_weekly.nc")
            # coarsen data
            data_2x2=data.coarsen(latitude=2, longitude=2).sum(['latitude','longitude'])
            data_2x2.to_netcdf(f"{dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_2x2_weekly.nc")
            print(f"saved {dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_2x2_weekly.nc")
            
            data_4x4=data.coarsen(latitude=4, longitude=4).sum(['latitude','longitude'])
            data_4x4.to_netcdf(f"{dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_4x4_weekly.nc")
            print(f"saved {dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_4x4_weekly.nc")
            del data, data_2x2, data_4x4

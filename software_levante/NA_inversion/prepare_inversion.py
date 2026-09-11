# use conda environment ('inversion' for Sarah Grandke)

#   processes Flexpart runs 
#       calculates last position of Flexpart particles
#       interpolated concentration for each particle --> background concentration for each measurement
#       interpolated concentration at measurement
#       remaining particles
#   calculates hourly FLEXPART footprints
#       weekly footprints
#       weekly footprints- scaled diurnal
#       weekly footprints additive diurnal
#   coarsens resolution of footprints

# To use it, adapt:
#   start_date, end_date
#   directory: parent directory for all footprints and evaluations
#   BG (background model) & directory to BG concentration fiels
#   directory to reference fluxes and prior fluxes
#   diurnal offset/scaling
#   region
#   resolution
#   FLEXPART directory
#   mesurement datasets

# individualize runs with
#   estimate_csv_folder_name, prep_footprints_folder
#   datasets used in inversion
#   selective funtions run
#       PREP_TM5_4DVAR_REF_FLUXES
#       CALCULATE_LAST_POSITIONS
#       PROCESS_FLEX_RUNS=True   
#       GET_HIGH_RES_FOOTPRINTS=True       
#       GET_TM5_4DVAR_SCALED_FOOTPRINTS
#       GET_WEEKLY_NO_SCALING=False 
#       GET_WEEKLY_ADDITIVE_DIURNAL_FOOTPRINTS=True    
#       COARSEN_HIGH_RES_FOOTPRINT=True  ,.... 

import xarray as xr
import os
import pandas as pd
import datetime as dt
import numpy as np
from dateutil.relativedelta import relativedelta
from utils import get_start_date_of_week, p
import warnings


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

def createLastPositionRegionDF(path, region):
    '''Get position of each released particle, when particle leaves region save as dataframe
    Args:
        path: path to Flexpart output directory
        region: [lat_min, lat_max, lon_min, lon_max]
    Returns: nothing, saves DF_region_last_positions.pkl'''
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
    mask=data.lon.where((data.lat>=region[0]) & (data.lat<=region[1]) &(data.lon>=region[2]) & (data.lon<=region[3])).notnull()
    # True where particle enters the region (False -> True)
    entry = mask & ~mask.shift(time=1, fill_value=False)
    # Keep only particles that ever have an entry
    has_entry = entry.any("time")
    if not has_entry.all():
        bad_particles = has_entry.particle.where(~has_entry, drop=True)
        print(f"WARNINGParticles with no entry into the region: {bad_particles.values}")
        warnings.warn(f"WARNINGParticles with no entry into the region: {bad_particles.values}", UserWarning)
    valid = data.lon.notnull()
    last_valid_time = (valid.isel(time=slice(None, None, -1)).idxmax("time"))
    # Last entry time for each particle
    first_time =(entry.isel(time=slice(None, None, -1)).idxmax("time"))
    first_time = first_time.where(has_entry,other=last_valid_time)

    data_df=data.where(data.time==first_time).to_dataframe().dropna() 

    # save df
    print(f'Saving DataFrame to: {path}DF_last_region_positions.pkl')
    data_df.reset_index().to_pickle(path+'DF_last_region_positions.pkl')



class measurement_dataset():
    """Measurement dataset for which Flexpart was run

    Args:
        name: name of dataset
        version: version of dataset
        flexpart_path: Flexpart output directory
        measurement_path: path to measurement directory  
        mesurement_filename: filename example 
        csv_dir: path where .csv files will be saved in subdirectories yyyy_mm
        num_parts: number of particles used per measurement in Flexpart
        measured: 'co2' (e.g. insitu), 'xco2' (e.g. GOSAT)
        height_levels: 'single point' (e.g. in-situ), 'multiple points' (e.g. TCCON), 'continuous profile' (e.g. GOSAT RemoTeC)
        pressure_weighted: measurement weighted with pressure profile \n \t (if column_measurement==True)
        averaging_kernel: measurement weighted with averaging kernel \n \t (if column_measurement==True)
    """
    
    def __init__(
            self, 
            name: str, 
            version: str, 
            flexpart_path: str, 
            measurement_dir: str, 
            measurement_filename_str: str,
            measurement_filename_type: str,
            csv_dir: str,
            num_parts: int, 
            measured: str,
            height_levels: str,
            pressure_weighted: bool = False, 
            averaging_kernel: bool = False,
            averaging_kernel_intergration_operator: bool =False
            ):
        
        self.name, self.version = name, version
        self.flexpart_path, self.measurement_dir, self.csv_dir = flexpart_path, measurement_dir, csv_dir
        self.measurement_filename_str=measurement_filename_str
        self.measurement_filename_type=measurement_filename_type
        self.num_parts=num_parts

        self.measured, self.height_levels  = measured, height_levels
        self.pressure_weighted=pressure_weighted
        self.averaging_kernel=averaging_kernel
        self.averaging_kernel_intergration_operator=averaging_kernel_intergration_operator


    def create_all_last_positions(self, start_date, end_date):
        '''
        Calculate last positions for every day and saves them in Release directory
        '''
        for d in pd.date_range(start_date,end_date):
            release_dir=f'{self.flexpart_path}/{d.strftime("%Y_%m")}/Release_{d.strftime("%Y%m%d")}/'
            measurement_path=f'{self.measurement_dir}/{d.strftime("%Y_%m")}/{self.measurement_filename_str}{d.strftime("%Y%m%d")}{self.measurement_filename_type}'
            # check if there are measurements for that day
            if os.path.isfile(measurement_path):
                # calculate last positions of particles
                last_positions_path=f"{release_dir}/DF_last_positions.pkl"
                # check if DF_last_positions exists, create if doesnt exist
                if not os.path.isfile(last_positions_path):
                    print(f'Current date: {d.strftime("%Y%m%d")}')
                    print('create DF_last_positions.pkl')
                    createLastPositionsDF(release_dir)
    
    def create_all_last_region_positions(self, start_date, end_date):
        '''
        Calculate last region positions for every day and saves them in Release directory
        '''
        for d in pd.date_range(start_date,end_date):
            release_dir=f'{self.flexpart_path}/{d.strftime("%Y_%m")}/Release_{d.strftime("%Y%m%d")}/'
            measurement_path=f'{self.measurement_dir}/{d.strftime("%Y_%m")}/{self.measurement_filename_str}{d.strftime("%Y%m%d")}{self.measurement_filename_type}'
            # check if there are measurements for that day
            if os.path.isfile(measurement_path):
                # calculate last positions of particles
                last_positions_path=f"{release_dir}/DF_last_region_positions.pkl"
                # check if DF_last_region_positions exists, create if doesnt exist
                if not os.path.isfile(last_positions_path):
                    print(f'Current date: {d.strftime("%Y%m%d")}')
                    print('create DF_last_region_positions.pkl')
                    createLastPositionRegionDF(release_dir,[12, 56, -134, -62])

    def process_flexpart_runs(self, start_date, end_date, BG,BG_molefrac_bg_str, BG_molefrac_dir):
        '''
        calculate last positions, calculate background and TM5-4DVar values for flexpart runs of datasets 
        Args:
            start_date, end_date: dt.date() defining the time period
            BG: name of background model
            BG_molefrac_bg_str: submodel/version of background model
            BG_molefrac_dir: path to background model molefraction concentrations
            
            optional:
plim:   upper boundary of particle releases in hPa, defaults to 100hPa
        Returns:
            Nothing, will create insitu/TM5-4DVar_estimate/yyyy_mm and RemoTeCv240/TM5-4DVar_estimate/yyyy_mm subdirectories to save .csv files
        '''
        print('calc last positions')
        self.create_all_last_positions(start_date,end_date)
        # calculate background concentrations
        print('calc background values and interpolation values')
        for d in pd.date_range(start_date,end_date):
            print(f'Current date: {d.strftime("%Y%m%d")}')
            release_dir=f'{self.flexpart_path}/{d.strftime("%Y_%m")}/Release_{d.strftime("%Y%m%d")}/'
            measurement_path=f'{self.measurement_dir}/{d.strftime("%Y_%m")}/{self.measurement_filename_str}{d.strftime("%Y%m%d")}{self.measurement_filename_type}'
            print(measurement_path)
            # check if there are measurements for that day
            if os.path.isfile(measurement_path):
                if not os.path.isdir(f"{self.csv_dir}/{d.strftime('%Y_%m')}"):
                    print(f"create directory {self.csv_dir}/{d.strftime('%Y_%m')}")
                    os.makedirs(f"{self.csv_dir}/{d.strftime('%Y_%m')}")
                s_csv_path=f"{self.csv_dir}/{d.strftime('%Y_%m')}/{self.measured}_bg_{d.strftime('%Y%m%d')}.csv"
                # calc_background
                self.calc_background(release_dir=release_dir, measurement_path=measurement_path, s_csv_path=s_csv_path, BG_molefrac_dir=BG_molefrac_dir,interp_method='linear', col_name=f'{BG}_{BG_molefrac_bg_str}_background')
                # calc interpolation to background molefraction data
                self.calc_interpolated_BG(measurement_path,s_csv_path, BG_molefrac_dir,date=d, BGstr=f'{BG}_{BG_molefrac_bg_str}')
            else:
                print(f'No measurements given for {d.strftime("%Y%m%d")}')

    def calc_background(self,release_dir, measurement_path, s_csv_path, BG_molefrac_dir, interp_method='linear',col_name=''):
        """Calculate co2 background for in-situ measurements from dataframe containing last positions of all particles, path to TM5 data with pressure at boundaries 
        Args:
            BG_molefrac_dir: path to directory where concentration data is saved, needs to include pressure at boundaries
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
            col_name=f'background_{interp_method}'
        # read last positions data
        last_positions_df=pd.read_pickle(f"{release_dir}DF_last_positions.pkl")
        # check if pressure exists in df
        if not 'prs' in last_positions_df.columns:
            # caluculate height from pressure
            last_positions_df['prs']=p(last_positions_df.z)*100 # p[hPa]->p[Pa]
        
        # for last positions file, get timerange, read respective TM5 files
        date_min=np.min(last_positions_df.time).date()
        date_max=np.max(last_positions_df.time).date()
        # read molefrac data for that time range
        bg_list=[]
        for date in pd.date_range(date_min,date_max):
            date_str=date.strftime('%Y%m%d')
            bg_data=xr.open_dataset(f'{BG_molefrac_dir}/xco2_mean_{date_str}.nc')
            bg_list.append(bg_data)
        bg_data=xr.concat(bg_list, dim='times', data_vars='all')

        # add boundaries as coordinate 
        bg_data['boundaries']=bg_data.boundaries
        bg_data=bg_data.squeeze()

        # xarray with only particle number as dimension
        data=last_positions_df.set_index('particle').to_xarray()
        # get value of CT data based on particle time, lat and lon using interp_method
        bg_interp=bg_data.interp(times=data.time, latitude=data.lat, longitude=data.lon,method=interp_method)
        # get level the particles are in based on presssure
        # idxmax gives index of first time condition is true, -1 to get level nuber from upper boundary of level
        p_levels=((bg_interp.p_boundary<data.prs).idxmax(dim='boundaries')-1).reset_coords(names=['times','latitude','longitude'])      
        # set values of -1 to 0 (-1 if pressure lower than lowest boundary)
        p_levels['boundaries']=p_levels.boundaries.where(p_levels.boundaries != -1, 0)        # returns value from p_levels[boundary where condition is True, otherwise fills in other value (here 1)]
        # get co2 of respective level
        background=bg_interp.mix.sel(levels=p_levels.boundaries)
        # separate different sounding positions
        background_part=background.coarsen(particle=self.num_parts).mean()
        bg_interp.close()
        # get co2 value for nearest time & position for each insitu measurement position
        # read measurement data
        if self.measured=='co2' and self.measurement_filename_type=='.csv':
            # read measurement data
            measurement_data=pd.read_csv(measurement_path)        # , index_col=0
            # add background values to dataframe
            
        elif self.measured=='xco2':
            if self.averaging_kernel:
                measurement_data=xr.open_dataset(measurement_path)
                measurement_data['xco2_retrieved']=measurement_data.xco2
                if self.averaging_kernel_intergration_operator:
                    Axa=(measurement_data.integration_operator*measurement_data.xco2_averaging_kernel* measurement_data.co2_profile_apriori).sum('layer_dim')
                else:
                    Axa=(measurement_data.xco2_averaging_kernel/measurement_data.xco2_averaging_kernel.sum('layer_dim')*len(measurement_data.layer_dim)* measurement_data.co2_profile_apriori).mean('layer_dim')
                measurement_data['xco2']=measurement_data.xco2-measurement_data.co2_column_apriori+Axa
                if self.height_levels=='continuous profile':
                    measurement_data=measurement_data.drop_dims(['layer_dim','level_dim'])
                elif self.height_levels=='multiple points':
                    measurement_data=measurement_data.drop_dims(['layer_dim'])
                measurement_data=(measurement_data.to_dataframe().reset_index()).drop(columns=measurement_data.dims.keys())
        measurement_data[col_name]=background_part.values
        measurement_data.to_csv(s_csv_path, index=None)
        return

    def calc_interpolated_BG(self,measurement_path, s_csv_path, BG_molefrac_dir,date,BGstr=''):
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
            val = bg_data.xco2.interp(
                times=row['time'],
                latitude=row['latitude'],
                longitude=row['longitude'],
                method='linear'  # or 'nearest' if needed
            )
            return val.values.item()  # extract scalar
        def interpolate_co2_point(row):
            temp = bg_data.interp(
                times=row['time'],
                latitude=row['latitude'],
                longitude=row['longitude'],
                method='linear'  # or 'nearest' if needed
            )
            # surface pressure
            psurf=temp.pressure.item()
            # get pressure at intake_height
            p_intake=p(row['intake_height[magl]'], p0=psurf)
            val = temp.mix.swap_dims({'levels':'p_level'}).sel(p_level=p_intake, method="nearest")
            return val.values.item()  # extract scalar
        
        def interpolate_co2_at_averaging_kernel_continuous_profile(data):
            Axco2=[]
            for i in range(0,len(data.xco2)):
                bg_temp=bg_data.sel(latitude=data.latitude[i],longitude=data.longitude[i],times=data.time[i], method='nearest')
                data_temp=data.isel(sounding_dim=i)
                data_temp['pressure_levels']=data_temp.pressure_levels*100 # convert hPa in Pa
                xco2=0
                marker=False
                for k in reversed(range(len(data_temp.layer_dim))):
                    #print(f'{k=}')
                    layer_top_pressure=data_temp.pressure_levels.values[k+1]
                    layer_bottom_pressure=data_temp.pressure_levels.values[k]
                    A_i=data_temp.xco2_averaging_kernel.values[k]
                    if max(bg_temp.p_boundary.values)<layer_bottom_pressure:
                        marker=True
                        lower_level=0
                    else:
                        lower_level=len(bg_temp.where(bg_temp.p_boundary>=layer_bottom_pressure,drop=True).p_boundary)-1
                    #print(f'{lower_level=}')
                    if max(bg_temp.p_boundary.values)<layer_top_pressure:
                        upper_level=0
                    else:
                        upper_level=len(bg_temp.where(bg_temp.p_boundary>=layer_top_pressure,drop=True).p_boundary)-1
                    for bg_level in range(lower_level,upper_level+1):
                        p_1=min(bg_temp.p_boundary.values[bg_level],layer_bottom_pressure)
                        p_2=max(bg_temp.p_boundary.values[bg_level+1],layer_top_pressure)
                        p_diff = p_1-p_2
                        co2=bg_temp.isel(levels=bg_level).mix.values
                        xco2+=co2*A_i*p_diff  
                if marker:
                    xco2=xco2/(max(bg_temp.p_boundary.values)-min(data_temp.pressure_levels.values))  
                else:
                    xco2=xco2/(max(data_temp.pressure_levels.values)-min(data_temp.pressure_levels.values))
                Axco2.append(xco2)
            return Axco2
        
        # read data
        bg_data=xr.open_dataset(f'{BG_molefrac_dir}/xco2_mean_{date.strftime("%Y%m%d")}.nc').squeeze()
        if os.path.isfile(s_csv_path):  # check if file exists, there are some days without measurements
            csv_data=pd.read_csv(s_csv_path, parse_dates=['time'])
            # Apply interpolation function to each row in the DataFrame
            if self.measured == 'co2':
                bg_data['p_level']=((bg_data.p_boundary.sel(boundaries=slice(0,max(bg_data.boundaries.values)))+bg_data.p_boundary.sel(boundaries=slice(1,max(bg_data.boundaries.values)+1)))/2/100).swap_dims({'boundaries':'levels'})
                bg_data=bg_data.assign_coords(p_level=bg_data.p_level)
                csv_data[BGstr+'_interpolated_co2'] = csv_data.apply(interpolate_co2_point, axis=1)

            elif self.measured == 'xco2':  
                csv_data[BGstr+'_interpolated_xco2'] = csv_data.apply(interpolate_xco2_point, axis=1)
                measurement_data=xr.open_dataset(measurement_path)
                if self.averaging_kernel:
                    if self.height_levels=='multiple points':
                        print('not yet written!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    elif self.height_levels=='continuous profile':
                        csv_data[BGstr+'_interpolated_A*co2']=interpolate_co2_at_averaging_kernel_continuous_profile(measurement_data)
                    else:
                        print('you need to write that yourself !!!!')
            # save csv files
            csv_data.to_csv(s_csv_path, index=None)
        # insitu
        else: 
            print(f'no file {s_csv_path}')

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
        data=xr.concat(data_list, dim='time', join='outer')
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
    TM5_end_date=dt.date(2023,1,8)
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

     #needed to add this in since a CategoricalDtype with dates as categories can't be read as a data type? (S.G. 24.10.25)
    new_time_index = pd.to_datetime(ds_mean.indexes['time'].astype(str))
    ds_mean = ds_mean.assign_coords(time=new_time_index)
    if hasattr(ds_mean['time'], 'encoding'):
        for key in list(ds_mean['time'].encoding.keys()):
            # remove any suspicious encoding entries related to pandas/categorical
            if key in ('dtype', 'categories', 'pandas_type', 'pandas_version'):
                del ds_mean['time'].encoding[key]

    print(f'saving data to {spath}')
    ds_mean.to_netcdf(spath)
    return

# main
if __name__ == "__main__":  
    start_date, end_date =dt.date(2020,10,1), dt.date(2022,3,31)    # time period of measurements
    # data paths
    # parent directory, in which all Flexpart runs, measurement information is collected
    directory='/work/bb1170/RUN/b383736/data/Flexpart_2021/'
    

    # Background model:
    BG='TM5'   # 'TM5' or 'CAMS' molefractions with selected dataset as subdirectory
    #BG='CAMS'
    # directory to molefraction concentrations
    if BG=='TM5':
        BG_molefrac_dir=directory+'/TM54DVar/TM5_molfractions/'
        BG_molefrac_bg_str = 'RemoTeC_2.4.1+IS-land_ocean_bc'   # specifying str of folder to use
        
    elif BG=='CAMS':
        BG_molefrac_dir='/work/bb1170/RUN/b383736/data/CAMS/2020-2022/'
        BG_molefrac_bg_str = 'satellite'        # specifying str of folder to use

    # reference fluxes of TM5-4DVar and prior flux
    TM5flux_dir=directory+'/TM54DVar/TM54DVar_fluxes/'
    TM5_str = 'RemoTeC_2.4.1+IS'         # TM5-4DVar dataset, chose from ['RemoTeC_2.4.0+IS', 'ACOS+IS','IS', 'prior']

    prior_diurnal_path=directory+'/TM54DVar/high_res_total_scaling_RemoTeC+IS_scaling.nc' # path to diurnal offset or diurnal scaling data
    region=[12,56,-134,-62]   
    resolution=[2]      # resolutions for which inversion will be run

    # FLEXPART
    flex_dir=directory+'/Flexpart/' # path to directory with subdirectories for datasets
    # measurements directory
    meas_dir=directory+'/measurements/'

    # csv folder in which to save background concentrations, estimates from the background concentration, ... 
    estimate_csv_folder_name='/TM5-4DVar_estimate_3years/'
    prep_footprints_folder='/prep_footprints_TM5_3_years/' 
  

    # measurement datasets:
    gosat=measurement_dataset(
        name = 'RemoTeC', 
        version = '2.4.1', 
        flexpart_path = flex_dir+'/'+'RemoTeCv240',
        measurement_dir = meas_dir+'/'+'RemoTeC',
        measurement_filename_str = 'RemoTeCv2.4.1_',
        measurement_filename_type = '.nc',
        csv_dir = flex_dir+'/'+'RemoTeCv240'+'/'+estimate_csv_folder_name,
        num_parts = 40000,
        measured = 'xco2',
        height_levels = 'continuous profile',
        pressure_weighted = True,
        averaging_kernel = True)
    
    insitu=measurement_dataset(
        name='insitu',
        version='GLOBALVIEWplus_v10.1',
        flexpart_path=flex_dir+'insitu',
        measurement_dir=meas_dir+'Obspack',
        measurement_filename_str = 'ISpositions_',
        measurement_filename_type = '.csv',
        csv_dir=flex_dir+'insitu'+estimate_csv_folder_name,
        num_parts=40000,
        measured = 'co2',
        height_levels = 'single point')
    
    tccon=measurement_dataset(
        name='TCCON',
        version='',
        flexpart_path=flex_dir+'TCCON',
        measurement_dir=meas_dir+'TCCON',
        measurement_filename_str = 'TCCON_',
        measurement_filename_type = '.nc',
        csv_dir=flex_dir+'TCCON'+estimate_csv_folder_name,
        num_parts=40000,
        measured = 'xco2',
        height_levels = 'multiple points',
        pressure_weighted=True,
        averaging_kernel=True,
        averaging_kernel_intergration_operator=True)

    datasets=[gosat,tccon]
    # select funtions to run
    PREP_TM5_4DVAR_REF_FLUXES=False # calculate weekly prior fluxes, cut them to inversion region and adjust resolution to the one needed in inversion
    CALCULATE_LAST_POSITIONS=False   # sometimes this needs more computational memory than available on shared node, so this can be run seprately from the processing on a compute node   
    PROCESS_FLEX_RUNS=False     # calculates last positions if needed, background concentrations, interpolated values from the background dataset, and remaining fraction of particles after the flexpart run
      

    GET_HIGH_RES_FOOTPRINTS=False        # combines footprint data for all days in a month with the measurement data
    GET_TM5_4DVAR_SCALED_FOOTPRINTS=False    # will also calculate unscaled weekly footprints, optionally: remove hourly footprint files for storage reasons
    GET_WEEKLY_NO_SCALING=False         # use if no scaling should be applied to the footprints, adapt path in COARSEN_HIGH_RES_FOOTPRINT in that case
    GET_WEEKLY_ADDITIVE_DIURNAL_FOOTPRINTS=False     # use if additive diurnal cycle should be substracted

    COARSEN_HIGH_RES_FOOTPRINT=True     #False     # if True, adapt coarsen_footprint_in_folder
    coarsen_footprint_in_folder='scaled_weekly'

    Delete_High_Res_Footprints = False  # delete hourly high resolution footprints to reduce storage space

    if PREP_TM5_4DVAR_REF_FLUXES:
        # cuts specified region, coarsens to desired resolution, saves weekly files
        # for dataset selected for background calculation and prior
        for temp in [TM5_str, 'prior']:
            flux_path=f'{TM5flux_dir}/flux_1x1_{temp}.nc'
            # cut fluxes to desired region for 1x1 and res resolution
            for res in [1]+resolution: 
                prep_TM5_4DVar_flux(flux_path, region, res)
                # get weekly fluxes for resxres 
                get_weekly_TM5_4DVarflux(f"{flux_path.replace('_1x1_', f'_{res}x{res}_')[:-3]}_cut.nc")
     
    if CALCULATE_LAST_POSITIONS:
        for dataset in datasets:
            print(dataset.name)
            #dataset.create_all_last_positions(start_date,end_date)
            dataset.create_all_last_region_positions(start_date,end_date)
    
    # calculate last positions (if not already calculated), calculate background concentrations and estimates/interpolations from background model
    if PROCESS_FLEX_RUNS:
        for dataset in datasets: # [insitu, gosat, tccon]:
            print(dataset.name)
            dataset.process_flexpart_runs(
                start_date=start_date,
                end_date=end_date,
                BG= BG,
                BG_molefrac_bg_str= BG_molefrac_bg_str,
                BG_molefrac_dir=BG_molefrac_dir+BG_molefrac_bg_str)
            get_frac_remaining_particles(
                dir_path=dataset.flexpart_path,
                csv_dir=dataset.csv_dir,
                csv_sdir=dataset.csv_dir, 
                file_str=dataset.measured+'_bg',
                start_date=start_date,
                end_date=end_date, 
                num_parts=dataset.num_parts)
            
    # get 1x1 hourly footprints with measurement data into one ds per month, needed for diurnal cacle
    if GET_HIGH_RES_FOOTPRINTS:
        for month_start in pd.date_range(start_date, end_date, freq='MS'):
            month_end= month_start+relativedelta(months=1, days=-1)
            print(month_start, month_end)
            dir_list=[f'/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/' for date in pd.date_range(month_start, month_end) if date in pd.date_range(month_start, month_end)]
            
            # for insitu data
            for dataset in datasets: # [insitu,gosat,tccon]:
                files=[dataset.flexpart_path+d+f for d in dir_list if os.path.isdir(dataset.flexpart_path+d) for f in os.listdir(dataset.flexpart_path+d) if f.startswith('grid_time_')]
                ds=[]
                for f in files:
                    temp=xr.open_dataset(f, decode_timedelta=True).sel(height=30)[['spec001_mr']].squeeze(dim='nageclass')  
                    temp['spec001_mr'] = (temp.spec001_mr /30*28.96/44*10**6).assign_attrs(description='Flexpart footprint in units necessary for inversion', units='ppm s m^2/kgCO2') # divide by layer height - s m^3/kg --> s m^2 /kg
                    temp['release_num']=(('pointspec'), temp.pointspec.values)
                    temp['release_day']=pd.to_datetime(f[-36:-28])
                    ds.append(temp)
                footprint=xr.concat(ds, dim='pointspec', data_vars='all', join='outer')
                footprint['pointspec']=footprint.pointspec.values
                data_list=[]
                for date in pd.date_range(month_start, month_end):    
                    # check if file exists for that date
                    path=f'{dataset.flexpart_path}/{estimate_csv_folder_name}/{date.strftime("%Y_%m")}/{dataset.measured}_bg_{date.strftime("%Y%m%d")}.csv'
                    if os.path.isfile(path):
                        data=pd.read_csv(path, parse_dates=['time'])
                        data_list.append(data)
                # add meas info to this dataset
                measurement_data=pd.concat(data_list)
                measurement_data=measurement_data.reset_index(names='release_num')
                # check that release dates match
                if (footprint.release_day.dt.date.values==measurement_data.time.dt.date).all():
                    # check that release numbers match 
                    if (footprint.release_num.values==measurement_data.release_num).all():
                        # rename columns
                        measurement_data.rename(columns={'latitude':'release_lat', 'longitude':'release_lon', 'time':'release_time'}, inplace=True)
                        cols=measurement_data.columns
                        for col in cols:
                            footprint[col]=(('pointspec'), measurement_data[col])
                        if not os.path.isdir(f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly'):
                            os.makedirs(f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly')
                        s_file_path=f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                        footprint.to_netcdf(s_file_path)
                        print(f'saved to {s_file_path}')
                    else:
                        print('Problem with Release numbers')
                else: 
                    print('Problem with Release days')

    # apply TM5-4DVar diurnal cycle scaling, add to weekly
    if GET_TM5_4DVAR_SCALED_FOOTPRINTS:
        print('reading scaling data')
        scaling_data=xr.open_dataset(prior_diurnal_path)
        for dataset in datasets:  #
            for month_start in pd.date_range(start_date, end_date, freq='MS'):
                month_end= month_start+relativedelta(months=1, days=-1)
                # read footprint data
                data_path=f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                print('reading footprint data')
                data=xr.open_dataset(data_path)
                # multiply with TM5-4DVar total scaling factors from prior
                data['spec001_mr_scaled']=(data['spec001_mr']*scaling_data.scaling_total).assign_attrs(description='Flexpart footprint in units necessary for inversion, scaled with total diurnal scaling factor from TM5-4DVar prior fluxes (sc_tot = bio_h/tot_mean + rest_mean/tot_mean)', units='ppm s m^2/kgCO2')
    
                # save dataset
                data_spath=f'{dataset.flexpart_path}/{prep_footprints_folder}/scaled_hourly/high_res_scaled_footprints_{month_start.strftime("%Y_%m")}.nc'
                if not os.path.isdir(f'{dataset.flexpart_path}/{prep_footprints_folder}/scaled_hourly'):
                    os.makedirs(f'{dataset.flexpart_path}/{prep_footprints_folder}/scaled_hourly')
                print(f'saving preprocessed hourly footprints to {data_spath}')
                data.to_netcdf(data_spath, mode="w")
                print('saving succesfull')

                # get weekly sum
                # Create the 7-day period bins
                period_bins = pd.date_range(start=get_start_date_of_week(month_start-dt.timedelta(days=10)), end=(get_start_date_of_week(month_end)+dt.timedelta(days=7)), freq="7D")
                # Cut the time array into the defined 7-day bins
                time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=pd.to_datetime(period_bins[:-1]))
                # Assign the new time bins as coordinates
                data=data.assign_coords(time=time_bins)
                # weekly sum
                # the timestamp corresponds to the first day of the week
                data=data.groupby("time").sum('time')
                data=data.assign_coords(time=data.time.astype("datetime64[ns]"))
                # save dataset
                spath=f"{dataset.flexpart_path}/{prep_footprints_folder}/scaled_weekly/high_res_scaled_footprints_{month_start.strftime('%Y_%m')}_weekly.nc"
                if not os.path.isdir(f'{dataset.flexpart_path}/{prep_footprints_folder}/scaled_weekly'):
                    os.makedirs(f'{dataset.flexpart_path}/{prep_footprints_folder}/scaled_weekly')
                print(f'saving preprocessed weekly footprints to {spath}')
                data.to_netcdf(spath, mode="w")
                print('saving succesfull')
                
                del data

                # optional, to limit storage space
                if Delete_High_Res_Footprints:
                    # remove old high_res hourly footprint file
                    print(f'deleting old hourly footprint file: {data_path}')
                    os.remove(data_path)
                
    if GET_WEEKLY_ADDITIVE_DIURNAL_FOOTPRINTS:
        diurnal_cycle_data=xr.open_dataset(prior_diurnal_path)
        for dataset in datasets: #[insitu,gosat,tccon]:   #
            for month_start in pd.date_range(start_date, end_date, freq='MS'):
                month_end= month_start+relativedelta(months=1, days=-1)
                # read footprint data
                data_path=f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                print('reading footprint data')
                data=xr.open_dataset(data_path)
                # multiply with TM5-4DVar total scaling factors from prior
                data['spec001_mr_diurnal_amplitude']=(data['spec001_mr']*diurnal_cycle_data.nee_amplitude*44/12*1e-3/(24*60*60)).assign_attrs(description='additive diurnal term from TM5-4DVar prior fluxes hourly_amplitude * Flexpart footprint', units='ppm s m^2/kgCO2')

                # get weekly sum
                # Create the 7-day period bins
                period_bins = pd.date_range(start=get_start_date_of_week(month_start-dt.timedelta(days=10)), end=(get_start_date_of_week(month_end)+dt.timedelta(days=7)), freq="7D")
                # Cut the time array into the defined 7-day bins
                time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=pd.to_datetime(period_bins[:-1]))
                # Assign the new time bins as coordinates
                data=data.assign_coords(time=time_bins)
                # weekly sum
                # the timestamp corresponds to the first day of the week
                data=data.groupby("time").sum('time')
                data=data.assign_coords(time=data.time.astype("datetime64[ns]"))

                # calculating offset of measurements due to diurnal cycle
                diurnal_offset=data['spec001_mr_diurnal_amplitude'].stack(grid_box=("time","latitude", "longitude")).squeeze()
                diurnal_offset=diurnal_offset.sum('grid_box')
                data[f'real_{dataset.measured}']=data[dataset.measured]
                data[dataset.measured] = data[dataset.measured]-diurnal_offset

                # save dataset
                spath=f'{dataset.flexpart_path}/{prep_footprints_folder}/diurnal_amplitude_weekly/high_res_scaled_footprints_{month_start.strftime("%Y_%m")}_weekly.nc'
                if not os.path.isdir(f'{dataset.flexpart_path}/{prep_footprints_folder}/diurnal_amplitude_weekly'):
                    os.makedirs(f'{dataset.flexpart_path}/{prep_footprints_folder}/diurnal_amplitude_weekly')
                print(f'saving preprocessed weekly additive diurnal term  to {spath}')
                data.to_netcdf(spath, mode="w")
                print('saving succesfull')
                
                del data

                if Delete_High_Res_Footprints:
                    # remove old high_res hourly footprint file to save storage space
                    print(f'deleting old hourly footprint file: {data_path}')
                    os.remove(data_path)

    if GET_WEEKLY_NO_SCALING:
        for dataset in datasets:  #
            for month_start in pd.date_range(start_date, end_date, freq='MS'):
                month_end= month_start+relativedelta(months=1, days=-1)
                # read footprint data
                data_path=f'{dataset.flexpart_path}/{prep_footprints_folder}/hourly/high_res_footprints_{month_start.strftime("%Y_%m")}.nc'
                print('reading footprint data')
                data=xr.open_dataset(data_path)

                # get weekly sum
                # Create the 7-day period bins
                period_bins = pd.date_range(start=get_start_date_of_week(month_start-dt.timedelta(days=10)), end=(get_start_date_of_week(month_end)+dt.timedelta(days=7)), freq="7D")
                # Cut the time array into the defined 7-day bins
                time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=pd.to_datetime(period_bins[:-1]))
                # Assign the new time bins as coordinates
                data=data.assign_coords(time=time_bins)
                # weekly sum
                # the timestamp corresponds to the first day of the week
                data=data.groupby("time").sum('time')
                data=data.assign_coords(time=data.time.astype("datetime64[ns]"))

                # save dataset
                spath=f"{dataset.flexpart_path}/{prep_footprints_folder}/weekly/high_res_footprints_{month_start.strftime('%Y_%m')}_weekly.nc"
                if not os.path.isdir(f'{dataset.flexpart_path}/{prep_footprints_folder}/weekly'):
                    os.makedirs(f'{dataset.flexpart_path}/{prep_footprints_folder}/weekly')
                print(f'saving preprocessed weekly footprints to {spath}')
                data.to_netcdf(spath, mode="w")
                print('saving succesfull')
                del data
                if Delete_High_Res_Footprints:
                    # remove old high_res hourly footprint file to save storage space
                    print(f'deleting old hourly footprint file: {data_path}')
                    os.remove(data_path)
    # combine into one ds for entire time period,  coarsen to 2x2 and 4x4 resolution
    if COARSEN_HIGH_RES_FOOTPRINT:
        for dataset in datasets:#,gosat,tccon]:
            print(f'read data for {dataset.name}')
            dir_path=f'{dataset.flexpart_path}/{prep_footprints_folder}/{coarsen_footprint_in_folder}/'   
            # Preprocess function: only select first timestep for specific cols
            def preprocess_get_first_timestep(ds):
                exclude_vars = {"spec001_mr",'spec001_mr_diurnal_amplitude','spec001_mr_scaled'}  
                # Drop unwanted variables early
                for var in ds.data_vars:
                    if var in exclude_vars:
                        continue
                    # if var in ds:
                    if "time" in ds[var].dims:
                        ds[var] = ds[var].isel(time=0, drop=False)  # keep time dimension (size 1)  # keeps time as length-1 dim
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
            
            data.to_netcdf(f"{dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_weekly.nc")
            print(f"saved {dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_weekly.nc")
            # coarsen data
            for res in resolution:
                data_resxres=data.coarsen(latitude=res, longitude=res).sum(['latitude','longitude'])
                data_resxres.to_netcdf(f"{dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_{res}x{res}_weekly.nc")
                print(f"saved {dir_path}/high_res_scaled_footprints_{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_{res}x{res}_weekly.nc")
                del data_resxres
            del data

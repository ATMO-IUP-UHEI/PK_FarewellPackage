from os import listdir, makedirs
from os.path import isfile, join, exists, isdir
import xarray as xr
import yaml
import click
import numpy as np
import datetime as dt
import pandas as pd

###################################################
# run this script e.g. with 
# python create_GOSATpositions.py --gosat_config_path configs/gosat_config_RemoTeC240.yaml
###################################################

@click.command()
@click.option("--gosat_config_path", required=True, help="path to gosat_config.yaml file", type=str)
@click.option('--outfile','-o',help='filename for output', type=str)
def main(gosat_config_path,outfile):
    GOSAT_dir,RemoTeC_version, Lat_min, Lat_max, Long_min, Long_max, outdir,startdate,enddate=read_gosat_config(gosat_config_path)
    if outdir is None:
        outdir=GOSAT_dir
        
    # RemoTeC v2.4.1
    if RemoTeC_version=='2.4.1':
        write_GOSATpositions_RemoTeC241(GOSAT_dir, Lat_min, Lat_max, Long_min, Long_max,outdir, startdate,enddate)
    
    # RemoTeC v2.4.0
    if RemoTeC_version=='2.4.0':
        start_date=dt.date(int(startdate[:4]),int(startdate[4:6]),int(startdate[6:]))
        end_date=dt.date(int(enddate[:4]),int(enddate[4:6]),int(enddate[6:]))
        write_GOSATpositions_RemoTeC240(GOSAT_dir, Lat_min, Lat_max, Long_min, Long_max,outdir,start_date,end_date)
    
def read_gosat_config(gosat_config_path):
    with open(gosat_config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        return (config['GOSAT_dir'],config['RemoTeC_version'], config['Lat_min'], config['Lat_max'], config['Long_min'], config['Long_max'],config['outdir'], config['startdate'],config['enddate'])

def write_GOSATpositions_RemoTeC241(GOSAT_dir, Lat_min, Lat_max, Long_min, Long_max,outdir, startdate,enddate, outfile='GOSATpositions'):  
    # GOSAT_dir: path to directory containing GOSAT files  
    # startdate,enddate: 
    GOSAT_files=[f for f in listdir(GOSAT_dir) if isfile(join(GOSAT_dir, f)) if f.endswith('.nc')]
    GOSAT_files.sort()
    if not startdate==None:
        GOSAT_files=[f for f in GOSAT_files if int(f[-11:-3]) >= int(startdate)]
    if not enddate==None:
        GOSAT_files=[f for f in GOSAT_files if int(f[-11:-3]) <= int(enddate)]
    for file in GOSAT_files:
        data=xr.open_dataset(GOSAT_dir+file)
        print(f'reading file: {GOSAT_dir}{file}')
        # only use data over land and nadir and cut box
        data=data.where((data.flag_landtype==0) & (data.flag_sunglint==0) & (data.latitude>Lat_min) & (data.latitude<Lat_max) & (data.longitude>Long_min) ,drop=True)
        # check if there is longitude < Long_min, if so cut and continue
        if np.min(data.longitude)<Long_max:
            data=data.where((data.longitude<Long_max), drop=True)
            # average if meas within 0.5° and within 5 minutes
            meas_id=1
            meas_id_array=np.zeros(len(data.sounding_dim))
            meas_id_array[0]=meas_id
            for i in range(1,len(data.sounding_dim)):
                dlon=np.abs(data.longitude.values[i]-data.longitude.values[i-1])
                dlat=np.abs(data.latitude.values[i]-data.latitude.values[i-1])
                dt=np.abs(np.datetime64(data.time.values[i])-np.datetime64(data.time.values[i-1]))
                if dlon>0.5 or dlat>0.5 or dt>np.timedelta64(5,'m'):
                    meas_id+=1
                meas_id_array[i]=meas_id        
            # print(meas_id_array)                
            # drop all variables with dependencies other than soundig_dim
            # TODO nachfragen variable gain
            data_df=data.drop_vars(['co2_profile_apriori','ch4_profile_apriori','dry_airmass_layer','pressure_levels','xco2_averaging_kernel','xch4_averaging_kernel','gain']).to_dataframe()
            data_df['meas_id']=meas_id_array
            data_df=data_df.groupby(['meas_id']).mean()
            date_str=pd.to_datetime(data.time.values[0]).strftime('%Y%m%d')
            filename=outdir+outfile+f'_{date_str}.csv' # add date string to filename
            data_df.to_csv(filename, columns=['time','latitude','longitude','xco2', 'xco2_err'],index=False)   #,mode='a',header=(not exists(filename)))  # append to file if already exists if meas for multiple day in one .csv file
            print(f'saved sounding positions for region in: {filename}')
    return

def write_GOSATpositions_RemoTeC240(GOSAT_dir, Lat_min, Lat_max, Long_min, Long_max,outdir,start_date,end_date,outfile='RemoTeCv2.4.0', AVERAGE_MEAS=False): 
    ''' 
    Args:
        GOSAT_dir: path to directory containing yearly GOSAT RemoTeC v2.4.0 files
        Lat_min, Lat_max, Long_min, Long_max: region for which GOSAT files should be used
        outdir: path to output directories with subdirectories '%Y_%m'
        start_date, end_date: start and end date as dt.date
        outfile: string of output filename_%Y%m%d, defaults to 'RemoTeCv2.4.0'
        AVERAGE_MEAS: set to True if want to average over measurements close in space&time, might need to adjust limits below, defaults to False to dave single measurements
    Returns:
        nothing, saves gosat data as .csv files
    '''
    # read files for relevant years
    if start_date.year == end_date.year:
        ds=xr.open_dataset(f'{GOSAT_dir}/fp_short_fil_corr_201202_{start_date.year}.nc')
    else:
        ds_list=[]
        for y in range(start_date.year, end_date.year+1):
            ds=xr.open_dataset(f'{GOSAT_dir}/fp_short_fil_corr_201202_{y}.nc')
            ds_list.append(ds)
        ds=xr.concat(ds_list, dim='sounding_dim')
    ds=ds.assign_coords(time=ds.time, latitude=ds.latitude, longitude=ds.longitude).drop_dims(['layer_dim', 'level_dim'])
    # only use data over land and nadir and cut box
    ds=ds.where((ds.flag_landtype==0) & (ds.flag_sunglint==0) & (ds.latitude>Lat_min) & (ds.latitude<Lat_max) & (ds.longitude>Long_min) & (ds.longitude<Long_max),drop=True)
    val=0
    tot_val=0
    for date in pd.date_range(start_date, end_date):
        ds_sel=ds.where(ds.time.dt.date==date.date(), drop=True)
        if ds_sel.sounding_dim.size>0:
            data_df=ds_sel[['time','latitude','longitude','xco2', 'xco2_err']].to_dataframe()
            # check if output dir exists
            if not isdir(f"{outdir}{date.strftime('%Y_%m')}"):
                makedirs(f"{outdir}{date.strftime('%Y_%m')}")
                print(f"created directory: {outdir}{date.strftime('%Y_%m')}")
            filename=f"{outdir}{date.strftime('%Y_%m')}/{outfile}_{date.strftime('%Y%m%d')}.csv" # add date string to filename
            
            # average over measurements that are close together
            if AVERAGE_MEAS:
                meas_id=1
                meas_id_array=np.zeros(len(ds_sel.sounding_dim))
                meas_id_array[0]=meas_id
                for i in range(1,ds_sel.sounding_dim.size):
                    dlon=np.abs(ds_sel.longitude.values[i]-ds_sel.longitude.values[i-1])
                    dlat=np.abs(ds_sel.latitude.values[i]-ds_sel.latitude.values[i-1])
                    dt_val=np.abs(np.datetime64(ds_sel.time.values[i])-np.datetime64(ds_sel.time.values[i-1]))
                    # average if measurements within 0.15° and within a minute of each other
                    if dlon>0.15 or dlat>0.15 or dt_val>np.timedelta64(1,'m'):
                        meas_id+=1
                    meas_id_array[i]=meas_id
                data_df['meas_id']=meas_id_array
                data_df=data_df.groupby(['meas_id']).mean()
                filename=f"{outdir}{outfile}_{date.strftime('%Y%m%d')}_mean.csv" # add date string to filename
            print(f'saved sounding positions for region to: {filename}')
            data_df.to_csv(filename, index=False)

#########################################################################################

if __name__ == '__main__':
    main()
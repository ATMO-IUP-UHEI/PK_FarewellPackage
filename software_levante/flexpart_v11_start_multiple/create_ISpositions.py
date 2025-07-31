import os
import xarray as xr
import yaml
import click
import numpy as np
import datetime as dt
import pandas as pd
###################################################
# 
# run this script e.g. with 
# python create_ISpositions.py --config_path configs/IS_config.yaml
###################################################

@click.command()
@click.option("--config_path", required=True, help="path to IS_config.yaml file", type=str)
def main(config_path):
    is_file, outdir,startdate,enddate=read_IS_config(config_path)
    start_date=dt.date(int(startdate[:4]),int(startdate[4:6]),int(startdate[6:]))
    end_date=dt.date(int(enddate[:4]),int(enddate[4:6]),int(enddate[6:]))
    write_ISpositions(is_file, start_date, end_date, outdir)
    
def read_IS_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        return (config['path'], config['outdir'], config['startdate'],config['enddate'])

def write_ISpositions(is_file, start_date, end_date, outdir, outfile='ISpositions'):
    ''' 
    Args:
        is_file: path to .nc file containing obspack insitu data
        start_date, end_date: start and end date as dt.date
        outdir: path to output directories with subdirectories '%Y_%m', creates those if dont exist
        outfile: string of output filename_%Y%m%d, defaults to 'ISpositions'
    Returns:
        nothing, saves gosat data as .csv files
    '''
    # read meas data from file
    obs_data=xr.open_dataset(is_file)
    # select date
    for date in pd.date_range(start_date, end_date):
        # print(date.date())
        # check if mesurements exist for that day
        if date.date() in obs_data.time.dt.date:
            # select date
            obs_sel=obs_data.where(obs_data.time.dt.date==date.date(), drop=True)
            # Get indexes of non-NaN values for the target variable obs_sel.value
            indexes = np.array(np.where(~np.isnan(obs_sel.value.values))).T
            
            # select times, files, location etc 
            # save to .csv file
            num_releases=len(indexes)
            times_sel=[pd.to_datetime(obs_sel.time[indexes[i][1]].item()) for i in range(0,num_releases)]
            files_sel=[obs_sel.file[indexes[i][0]].item() for i in range(0,num_releases)]
            lat=[obs_sel.latitude.sel(time=times_sel[i], file=files_sel[i]).item() for i in range(0,num_releases)]
            lon=[obs_sel.longitude.sel(time=times_sel[i], file=files_sel[i]).item() for i in range(0,num_releases)]
            co2_val=[obs_sel.value.sel(time=times_sel[i], file=files_sel[i]).item()*1e6 for i in range(0,num_releases)] # mol/mol -> ppm molefraction
            elevation=[obs_sel.elevation.sel(time=times_sel[i], file=files_sel[i]).item() for i in range(0,num_releases)]
            intake_height=[obs_sel.intake_height.sel(time=times_sel[i], file=files_sel[i]).item() for i in range(0,num_releases)]

            df=pd.DataFrame(np.array([times_sel, files_sel, lat, lon, co2_val, elevation, intake_height]).T, columns=['time', 'file', 'latitude', 'longitude', 'co2_val[ppm]', 'elevation[masl]', 'intake_height[magl]'])
            # check if '%Y_%m' subdirectory exists
            if not os.path.isdir(f'{outdir}/{date.strftime("%Y_%m")}'):
                os.mkdir(f'{outdir}/{date.strftime("%Y_%m")}')
            spath=f'{outdir}/{date.strftime("%Y_%m")}/{outfile}_{date.strftime("%Y%m%d")}.csv'
            df.to_csv(spath)
            print(f'saved to {spath}')
    return

#########################################################################################

if __name__ == '__main__':
    main()
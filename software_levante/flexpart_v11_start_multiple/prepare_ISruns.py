import os
import xarray as xr
import numpy as np
from distutils.dir_util import copy_tree
import yaml
import create_part_init
import click
import pandas as pd
import datetime as dt
from os import listdir
from os.path import isfile, join
from prepare_GOSATruns import read_paths, read_sounding_positions, prepare_command, prepare_pathnames

###################################################
# run this script e.g. with 
# python prepare_ISruns.py --config_path configs/options_config_IS.yaml
###################################################

@click.command()
@click.option("--config_path", required=True, help="path to options_config.yaml file for IS releases",type=str)

def main(config_path):
    config=read_paths(config_path)
    options_dummy_path,output_path,input_path,available_path,sim_length,part_init_config_path = config['options_dummy_path'], config['output_path'], config['input_paths'], config['available_paths'],config['sim_length'], config['part_init_config_path']
    IS_positions_dir=config['IS_positions_dir']
    start_date=pd.to_datetime(config['startdate'], format='%Y%m%d').date()
    end_date=pd.to_datetime(config['enddate'], format='%Y%m%d').date()
    prepare_IS_releases(start_date, end_date, IS_positions_dir, output_path, options_dummy_path, input_path,available_path,sim_length,part_init_config_path)

def prepare_IS_releases(start_date, end_date, IS_positions_dir, output_path, options_dummy_path, input_path,available_path,sim_length,part_init_config_path):
    ''' 
    Args:
        start_date, end_date (dt.date): time range for releases
        IS_positions_dir: path to directory with the .csv files
        output_path: path to directory the Release directories should be saved into
        options_dummy_path:
    Returns:
        nothing, creates Release directories with initial conditions files, part_ic.nc
    '''
    # read Release specifications from part_init_config file
    part_init_config=create_part_init.read_config(part_init_config_path)
    num_part,nspecies,species_id,species_mass,kindz,num_layers = part_init_config['num_part'],part_init_config['nspecies'],part_init_config['species_id'],part_init_config['species_mass'],part_init_config['kindz'],part_init_config['num_layers']
    
    for date in pd.date_range(start_date, end_date):
        csv_path=f'{IS_positions_dir}/{date.strftime("%Y_%m")}/ISpositions_{date.strftime("%Y%m%d")}.csv'
        # check if mesurements exist for that day
        if os.path.isfile(csv_path):
            # read data
            df=pd.read_csv(csv_path, index_col=0, parse_dates=['time'], date_format='mixed')
            num_releases=df.index.size
            # start, end of simulation
            sim_max=np.max(df.time)+dt.timedelta(hours=4) # 
            sim_min=np.min(df.time)-dt.timedelta(days=sim_length) # 
            # particles released during 4h period, equaly spaced during that time
            dt_list=np.linspace(0,4*60*60, num=num_part)
            release_time_end=(sim_max-df.time).dt.total_seconds()     # distance of df.time to end of simulation in seconds
            
            if not (num_part % num_layers)==0:
                print('ERROR: number of particles can not be evenly distributed on number of layers')
            # adjust numpart so that number of particles per layer is equal
            part_per_layer=int(num_part/num_layers)
            # get particle ids, ranges from 0 to number of particles per measurement release
            particle_id=np.arange(0,num_part*num_releases,step=1)

            # release ids for one measurement
            rmeas=np.zeros(num_part)
            for j in range(0,num_layers):
                rmeas[j*part_per_layer:(j+1)*part_per_layer]=j+1   # separate into num_layers number of layers

            # create array of needed length
            # longitude 	Initial longitude of each particle 	1D-array of reals with dimension particle
            longitude = np.ones(num_part*num_releases)
            # latitude 	Initial latitude of each particle 	1D-array of reals with dimension particle
            latitude = np.ones(num_part*num_releases)
            # time 	Release time of each particle seconds after simulation start (IBDATE/IBTIME for forward runs, IEDATE/IETIME for backward runs, set in COMMAND) 	1D-array of integers with dimension particle
            time = np.ones(num_part*num_releases)
            # height 	Initial height of each particle (meter above reference level) 	1D-array of reals with dimension particle
            height = np.zeros(num_part*num_releases)
            # release 	Release ID of each particle, giving separate concentration fields for each ID when IOUTPUTFOREACHRELEASE in COMMAND is set 	1D-array of integers with dimension particle
            release=np.zeros(num_part*num_releases)

            # set values for each release
            for i in range(0,num_releases):
                longitude[i*num_part:(i+1)*num_part]=df.longitude[i]
                latitude[i*num_part:(i+1)*num_part]=df.latitude[i]
                time[i*num_part:(i+1)*num_part]=release_time_end[i]-dt_list
                height[i*num_part:(i+1)*num_part]=df['intake_height[magl]'][i]
                release[i*num_part:(i+1)*num_part]=rmeas+i*num_layers
            # same mass for all releases
            # mass 	Initial mass of each particle (kg) 	2D-array of reals with dimension species and particle
            mass=np.ones((1,num_part*num_releases))
            # mass[0] = species   # species
            mass[0] = species_mass/6.022*10**-26           # mass in kg; molarweight[g/mol]/N_A/10**3[g->kg]

            # create xarray
            ds=xr.Dataset(
                data_vars=dict(
                    longitude=(['particle'],longitude),
                    latitude=(['particle'],latitude),
                    height=(['particle'],height),
                    time=(['particle'],time),
                    mass=(['species','particle'],mass),
                    release=(['particle'],release)
                    ),
                    coords=dict(
                        particle=('particle',particle_id)
                    ),
                    attrs=dict(nspecies=(nspecies),
                            species=(species_id),
                            kindz=(kindz),
                            num_layers=(num_layers,'number of layers per release'),
                            num_releases=(num_releases,'number of IS releases'))
            )
            # create Release_directory
            release_outdir_day=f'{output_path}/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/'
            os.makedirs(release_outdir_day)
            # save dataset
            ds.to_netcdf(release_outdir_day+'part_ic.nc')
            print(f'Written file {release_outdir_day}part_ic.nc')
            
            # make options directory
            options_outpath_temp=f'{output_path}/{date.strftime("%Y_%m")}/config/options_{date.strftime("%Y%m%d")}/'
            os.makedirs(options_outpath_temp,exist_ok = True)
            # copy all files from options dummy directory
            copy_tree(options_dummy_path, options_outpath_temp)
            # write command file
            prepare_command(options_dummy_path,options_outpath_temp, sim_max.to_numpy().astype('datetime64[s]'), sim_min.to_numpy().astype('datetime64[s]'))
            # write pathnames file
            os.makedirs(f'{output_path}{date.strftime("%Y_%m")}/config/pathnames/',exist_ok = True)
            # pathnames_file=f'{output_path}config/pathnames/pathnames_{release_id}'    # pathnames files numbered
            pathnames_file=f'{output_path}{date.strftime("%Y_%m")}/config/pathnames/pathnames_{date.strftime("%Y%m%d")}'      # uses dateinstead of number
            prepare_pathnames(pathnames_file,options_outpath_temp, release_outdir_day,input_path,available_path)
    return

#########################################################################################

if __name__ == '__main__':
    main()
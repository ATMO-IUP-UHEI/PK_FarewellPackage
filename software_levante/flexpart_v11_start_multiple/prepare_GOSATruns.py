import os
import numpy as np
from distutils.dir_util import copy_tree
import yaml
import create_part_init
import click
import pandas as pd
import datetime as dt
from os import listdir
from os.path import isfile, join

###################################################
# run this script e.g. with 
# python prepare_GOSATruns.py --config_path configs/options_config_RemoTeC240.yaml
###################################################

@click.command()
@click.option("--config_path", required=True, help="path to options_config.yaml file",type=str)
def main(config_path):
    config=read_paths(config_path)
    options_dummy_path,output_path,input_path,available_path,sim_length,part_init_config_path= config['options_dummy_path'], config['output_path'], config['input_paths'], config['available_paths'],config['sim_length'],config['part_init_config_path']
    sounding_dir= config['satellites_position_path']
    start_date=pd.to_datetime(config['startdate'], format='%Y%m%d').date()
    end_date=pd.to_datetime(config['enddate'], format='%Y%m%d').date()
    # if no sounding path given, use release_ids and locations specified in options_config.yaml
    if sounding_dir==None:
        print('Using release_ids, times and locations specified in options_config.yaml')
        sim_start,release_ids, release_time=config['sim_start'], config['release_ids'], config['release_time']
        latitude, longitude=config['latitude'], config['longitude']
        num_releases=config['num_releases']
        sim_start_last=sim_start
        # check if num_releases is integer, if so create list of 
        if isinstance(num_releases, int):
            num_releases=[num_releases]*len(release_ids)
        # check that lat, lon, sim_start and release_time have right format
        for i in range(0,len(release_ids)):
            if not (len(latitude[i])==num_releases[i]):
                print('ERROR: wrong format for latitude')
                return
            if not (len(longitude[i])==num_releases[i]):
                print('ERROR: wrong format for longitude')
                return
            if not (len(sim_start[i])==2):
                print('ERROR: wrong format for sim_start')
                return
            if not (len(release_time[i])==num_releases[i]):
                print('ERROR: wrong format for release_time')
                return
    # otherwise read from satelight_sounding.csv
    else:   # with satelight_sounding path
        # creating empty lists for necessary variables
        sim_start, sim_start_last, release_time,latitude, longitude,num_releases =[],[],[],[],[],[]
        
        for date in pd.date_range(start_date, end_date):
        # for j in range(0,len(sounding_files)):
            soundings_path=f'{sounding_dir}/{date.strftime("%Y_%m")}/RemoTeCv2.4.0_{date.strftime("%Y%m%d")}.csv'
            # check that file exists
            if os.path.isfile(soundings_path):
                # with satelight_sounding path
                print(f'reading sounding positions from {soundings_path}')
                sounding_datetime,latitude_temp,longitude_temp=read_sounding_positions(soundings_path)
                latitude.append(latitude_temp)
                longitude.append(longitude_temp)
                num_releases.append(len(latitude_temp))
                # sounding_date,sounding_time,latitude,longitude=read_sounding_positions(soundings_path)
                
                # get min and max of sim_start to determine start and stop of simulation that covers sim_length for all releases
                sim_start_min=np.min(sounding_datetime).ceil('h')
                sim_start_max=np.max(sounding_datetime).ceil('h')
                
                # start of simulation
                sim_start.append([dt.datetime.strftime(sim_start_max, '%Y%m%d'),dt.datetime.strftime(sim_start_max, '%H%M%S')])
                # start for last/earliest release
                sim_start_last.append([dt.datetime.strftime(sim_start_min, '%Y%m%d'),dt.datetime.strftime(sim_start_min, '%H%M%S')])

                release_time.append(-np.floor((pd.to_datetime(sim_start_max)-sounding_datetime).dt.total_seconds()).astype(int))

    for i in range(0,len(release_time)):
        # get start and end of simulation as datetime64
        start=np.datetime64(f"{sim_start[i][0][:4]}-{sim_start[i][0][4:6]}-{sim_start[i][0][6:8]}T{sim_start[i][1][:2]}:{sim_start[i][1][2:4]}:{sim_start[i][1][4:6]}")
        stop=np.datetime64(f"{sim_start_last[i][0][:4]}-{sim_start_last[i][0][4:6]}-{sim_start_last[i][0][6:8]}T{sim_start_last[i][1][:2]}:{sim_start_last[i][1][2:4]}:{sim_start_last[i][1][4:6]}")+np.timedelta64(-sim_length, 'D')
        # print(start,stop)
        # get release time, lat,lon
        r_time=release_time[i]
        # print(r_time)
        lat=latitude[i]
        lon=longitude[i]
        num_r=num_releases[i]
        
        # make output directory
        # out_dir_temp=f'{output_path}Release_{release_id}/'  # use if Release directories should be numbered
        out_dir_temp=f'{output_path}/{sim_start[i][0][:4]}_{sim_start[i][0][4:6]}/Release_{sim_start[i][0]}/' 
        os.makedirs(out_dir_temp)
        # make pathnames directory
        for j in range(0,config['num_pathnames_dir']):
            pathnames_dir_temp=f'{output_path}/{sim_start[i][0][:4]}_{sim_start[i][0][4:6]}/config/pathnames_{j}'
            if not os.path.isdir(pathnames_dir_temp):
                os.makedirs(pathnames_dir_temp)

        # create part_init.nc file
        part_init_config=create_part_init.read_config(part_init_config_path)
        num_part, pmin,pmax,nspecies,species_id,species_mass,kindz,num_layers = part_init_config['num_part'],part_init_config['pmin'],part_init_config['pmax'],part_init_config['nspecies'],part_init_config['species_id'],part_init_config['species_mass'],part_init_config['kindz'],part_init_config['num_layers']
        create_part_init.create_part_init(out_dir_temp,num_part,num_layers,num_r, lat,lon,r_time,pmin,pmax,nspecies,species_id,species_mass,kindz)

        # make options directory
        # options_outpath_temp=f'{output_path}config/options_{release_id}/'
        options_outpath_temp=f'{output_path}/{sim_start[i][0][:4]}_{sim_start[i][0][4:6]}/config/options_{sim_start[i][0]}/'
        os.makedirs(options_outpath_temp,exist_ok = True)
        # copy all files from options dummy directory
        copy_tree(options_dummy_path, options_outpath_temp)
        # write command file
        prepare_command(options_dummy_path,options_outpath_temp, start, stop)
        # write pathnames file
        # pathnames_file=f'{output_path}config/pathnames/pathnames_{release_id}'    # pathnames files numbered
        pathnames_file=f'{output_path}/{sim_start[i][0][:4]}_{sim_start[i][0][4:6]}/config/pathnames_{i%config["num_pathnames_dir"]}/pathnames_{sim_start[i][0]}'      # uses dateinstead of number
        prepare_pathnames(pathnames_file,options_outpath_temp, out_dir_temp,input_path,available_path)
        
def read_sounding_positions(csv_path):
    # reads sounding time and positions from file
    # returns string of time as datetime, lat and lon
    position=pd.read_csv(csv_path, parse_dates=['time'], date_format='mixed')
    return position.time,position.latitude,position.longitude

def read_paths(config_path):
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        return config

# modified functions from /work/bb1170/RUN/b382762/software/ColumnFLEXPART/columnflexpart/scripts/flexpart_safe_start.py
def prepare_command(options_dummy_path: str,options_outpath: str, start: np.datetime64, stop: np.datetime64):
    """Sets start and stop for simulation in COMMAND file in desired options directory

    Args:
        options_path (str): Absolute path of options directory
        start (str): Value for IEDATE and IETIME (the other way around due to backwards run)
        stop (str): Value for IBDATE and IBTIME (the other way around due to backwards run)
    """    
    start_date, start_time = datetime64_to_yyyymmdd_and_hhmmss(start)
    stop_date, stop_time = datetime64_to_yyyymmdd_and_hhmmss(stop)
    with open(os.path.join(options_dummy_path, "COMMAND"), "r") as f:
        com_file = ""
        for line in f:
            addition = line
            if line.startswith(' IEDATE'):
                addition = f" IEDATE={start_date},\n"
            if line.startswith(' IETIME'):
                addition = f" IETIME={start_time},\n"
            if line.startswith(' IBDATE'):
                addition = f" IBDATE={stop_date},\n"
            if line.startswith(' IBTIME'):
                addition = f" IBTIME={stop_time},\n"
            com_file = com_file + addition
    # print(com_file)
    with open(os.path.join(options_outpath, "COMMAND"), "w") as f:
        f.write(com_file)
def prepare_pathnames(pathnames_file: str, options_path: str, output_path: str, input_paths: list[str], available_files: list[str]):
    """Writes paths into pathnames file

    Args:
        pathnames_file (str): Absolute path of pathnames file to write into
        options_path (str): Absolute path of options directory
        output_path (str): Absolute path of output directory
        input_paths (str): Absolute path of input files
        available_files (str): Absolute path of available files
    """    
    content = ""
    content += f"{options_path}\n"
    content += f"{output_path}\n"
    content += f"{input_paths}\n"
    content += f"{available_files}\n"    
    with open(pathnames_file, "w") as f:
        f.write(content)

# copied from /work/bb1170/RUN/b382762/software/ColumnFLEXPART/columnflexpart/scripts/flexpart_safe_start.py
def set_start(start: np.datetime64, step: int, shift: int) -> np.datetime64:
    """Finds siutable value for start of simulation (backwards in time)

    Args:
        start (np.datetime64): start of release
        step (int): Allowed step sizes starting from shift (eg with 3: start is possible each 3rd hour)
        shift (int): Shift to start steps (if 1 and step 3 staring values are: 1,4,7,10...) 

    Returns:
        np.datetime64: starting value suitable with step and shift
    """    
    assert np.timedelta64(1, "D") % np.timedelta64(step, "h") == 0, f"24 hours have to be devisable by value of step. (Your value: {step})"
    date = start.astype("datetime64[D]")
    shifts_late = np.arange(np.timedelta64(shift, "h"), np.timedelta64(shift+25, "h"), np.timedelta64(step, 'h'), dtype="timedelta64[h]")
    shifts_early = np.arange(np.timedelta64(shift, "h"), np.timedelta64(shift-25, "h"), - np.timedelta64(step, 'h'), dtype="timedelta64[h]")[::-1]
    shifts = np.concatenate([shifts_early, shifts_late])
    start_vals = date + shifts
    ret = start_vals[start_vals > start][0].astype("datetime64[s]")
    return ret
def get_start_stop(releases_file: str, step: int, shift: int, sim_lenght: int) -> tuple[np.datetime64, np.datetime64]:
    """Reads out start and stop out of releases file and calculates suitable start and stop values for simulation.

    Args:
        releases_file (str): Path to releases file
        step (int): Allowed step sizes starting from shift (eg with 3: start is possible each 3rd hour)
        shift (int): Shift to start steps (if 1 and step 3 staring values are: 1,4,7,10...)
        sim_lenght (int): Simulation length in days

    Returns:
        tuple[np.datetime64, np.datetime64]: Start, stop of simulation (from point of view of simulation start later then stop)
    """    
    with open(releases_file) as f:
        lines = f.readlines()
    start_dates = [yyyymmdd_to_datetime64(line.split("=")[1].split(",")[0].replace(" ", ""))
        for line in lines if "IDATE2" in line]
    start_times = [hhmmss_to_timedelta64(line.split("=")[1].split(",")[0].replace(" ", "")) 
        for line in lines if "ITIME2" in line]
    stop_dates = [yyyymmdd_to_datetime64(line.split("=")[1].split(",")[0].replace(" ", ""))
        for line in lines if "IDATE1" in line]
    stop_times = [hhmmss_to_timedelta64(line.split("=")[1].split(",")[0].replace(" ", "")) 
        for line in lines if "ITIME1" in line]
    
    starts = np.array(start_dates) + np.array(start_times)
    stops = np.array(stop_dates) + np.array(stop_times)
    # max/min/- instead of min/max/+ since runs go backwards
    start = max(starts)
    stop = min(stops)
    start = set_start(start, step, shift)
    stop = stop - np.timedelta64(sim_lenght, "D")
    return start, stop

# copied from /work/bb1170/RUN/b382762/software/ColumnFLEXPART/columnflexpart/utils/utils.py
def yyyymmdd_to_datetime64(date_string: str) -> np.datetime64:
    """Convert string of form yyyymmdd to np.datetime64.

    Args:
        date_string (str): String of date to convert 

    Returns:
        np.datetime64: Converted date
    """    
    date = np.datetime64(f"{date_string[:4]}-{date_string[4:6]}-{date_string[6:]}")
    return date
def datetime64_to_yyyymmdd_and_hhmmss(time: np.datetime64) -> tuple[str, str]:
    """Convert np.datetime64 to strings of type yyyymmdd and hhmmss 

    Args:
        time (np.datetime64): time to convert

    Returns:
        tuple[str, str]: String for date and time
    """    
    string = str(time)
    string = string.replace("-", "").replace(":", "")
    date, time = string.split("T")
    return date, time
    
#########################################################################################

if __name__ == '__main__':
    main()
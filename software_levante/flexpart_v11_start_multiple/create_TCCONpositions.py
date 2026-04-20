import json
import datetime as dt
import os
import xarray as xr
import yaml
import pandas as pd
import numpy as np
from timezonefinder import TimezoneFinder
import click

###################################################
# run this script e.g. with 
# python create_TCCONpositions.py --TCCON_config_path configs/TCCON_config.yaml
###################################################

@click.command()
@click.option("--TCCON_config_path", required=True, help="path to TCCON_config.yaml file", type=str)
@click.option('--outfile','-o',help='filename for output', type=str)
def main(tccon_config_path,outfile):
    data_dir,sites_path, Lat_min, Lat_max, Long_min, Long_max, outdir,startdate,enddate=read_config(tccon_config_path)
    if outdir is None:
        outdir=dir
    start_date=dt.date(int(startdate[:4]),int(startdate[4:6]),int(startdate[6:]))
    end_date=dt.date(int(enddate[:4]),int(enddate[4:6]),int(enddate[6:]))
    sites=find_sites(data_dir,sites_path, Lat_min, Lat_max, Long_min, Long_max,outdir, start_date,end_date)
    data=reading_data(data_dir,sites,start_date,end_date)
    for day in pd.date_range(start_date,end_date):
        if day.date() in data.time.dt.date.values:
            print(day)
            data_sel=data.where(data.time.dt.date==day.date(),drop=True)
            date_str=pd.to_datetime(data_sel.time.values[0]).strftime('%Y%m%d')
            directory=f"{outdir}/{pd.to_datetime(data_sel.time.values[0]).strftime('%Y_%m')}/"
            if not os.path.exists(directory):
                os.makedirs(directory)
            filename=f"{directory}/TCCON_{date_str}.nc"
            # add date string to filename
            print(f'saved measurements for {day} to: {filename}')
            data_sel.to_netcdf(filename)

def read_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        return (config['data_dir'],config['sites_path'], config['Lat_min'], config['Lat_max'], config['Long_min'], config['Long_max'],config['outdir'], config['startdate'],config['enddate'])

def find_sites(dir,sites_path, Lat_min, Lat_max, Long_min, Long_max,outdir, startdate,enddate):
    with open(sites_path, 'r', encoding='utf-8') as file:
        data = json.load(file)
    print('selecting measurement sites')
    sites=[]
    format_string = "%Y-%m-%d"
    for i in range(len(data)):
        if (data[i]['latitude']>Lat_min) & (data[i]['latitude']<Lat_max):
            if (data[i]['longitude']>Long_min) & (data[i]['longitude']<Long_max):
                if (data[i]['end_date'] == None) or (startdate< dt.datetime.strptime(data[7]['end_date'], format_string).date()):
                    sites.append(data[i])
    return sites

def prep_data(data_dir,filename,start_date,end_date):
    print('prepping data')
    data=xr.open_dataset(data_dir+filename, decode_timedelta=True)
    data=data.where(((data.time.dt.date>=start_date)& (data.time.dt.date<=end_date)), drop=True)
    data['latitude']=data.lat
    data['longitude']=data.long
    data=data[['time','latitude','longitude','xco2', 'xco2_error','ak_xco2','ak_pressure','ak_altitude']]

    start_time=dt.time(12,)
    end_time=dt.time(16,)

    tf = TimezoneFinder()
    tz_name=tf.timezone_at(lat=data.latitude.values[0], lng=data.longitude.values[0])
    local_time = (data.indexes["time"].tz_localize("UTC").tz_convert(tz_name).tz_localize(None))
    data["local_time"] = ("time", local_time.to_numpy())
    data_sel=data.where((data.local_time.dt.time>=start_time)& (data.local_time.dt.time<= end_time),drop=True)
    group_key = data_sel.local_time.dt.floor("D")
    data_mean = data_sel.groupby(group_key).mean()

    local_fixed = group_key.values + np.timedelta64(start_time.hour, "h")

    # convert date → datetime64 and add hours
    local_fixed = (data_mean["floor"].values.astype("datetime64[ns]")+ np.timedelta64(start_time.hour, "h"))
    utc_time = (pd.DatetimeIndex(local_fixed).tz_localize(tz_name).tz_convert("UTC").tz_localize(None))
    data_mean = (data_mean.assign_coords(time=("floor", utc_time)).swap_dims({"floor": "time"}).drop_vars("floor"))
    return data_mean

def reading_data(data_dir,sites,start_date, end_date):
    firstFile=True
    for i, site in enumerate(sites):
        print(site)
        for file in os.listdir(data_dir):
            if file.startswith(site['site_id']):
                print(f'Reading file {file}')
                if firstFile:
                    data=prep_data(data_dir=data_dir,filename=file,start_date=start_date,end_date=end_date)
                    firstFile=False
                else:
                    temp=prep_data(data_dir=data_dir,filename=file,start_date=start_date,end_date=end_date)
                    data=xr.concat([data,temp], dim='time',data_vars='all')
                break
    return data


#########################################################################################

if __name__ == '__main__':
    main()
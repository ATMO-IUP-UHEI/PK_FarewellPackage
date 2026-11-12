# downloaded data from server with  sftp -P 1022 sftp0041-1-r@ftp.prd.dip.wur.nl
# Password: akU8%S2Fewr!PspDVBqd&s4Y9L&@C9wCd*B
# go into folder which files you want to download (cd,ls,...)
# download with :'get filename'
# exit with 'quit'
import xarray as xr
import datetime as dt
import pandas as pd
from utils import get_start_date_of_week, p
import numpy as np

start_date, end_date =dt.date(2009,1,1), dt.date(2022,12,31)  
analy_region_lat=[2,66]
analy_region_lon=[-144,-52]
resolution=int(2/0.25)
vars=['CO2'] # ['DM', 'C', 'CO2', 'CO', 'CH4', 'NMOC_g', 'H2', 'NOx', 'N2O', 'PM2.5', 'TPC', 'OC', 'BC', 'SO2', 'C2H6', 'CH3OH', 'C2H5OH', 'C3H8', 'C2H2', 'C2H4', 'C3H6', 'C5H8', 'C10H16', 'C7H8', 'C6H6', 'C8H10', 'Toluene_Lump', 'Higher_Alkenes', 'Higher_Alkanes', 'CH2O', 'C2H4O', 'C3H6O', 'NH3', 'C2H6S', 'HCN', 'HCOOH', 'CH3COOH', 'MEK', 'CH3COCHO', 'HOCH2CHO']
version='GFED5'

if version=='GFED5':
    GFED_path='/work/bb1170/RUN/b383736/data/GFED/GFED5.1/'
    GFED_NRT_path='/work/bb1170/RUN/b383736/data/GFED/GFED5NRT/'
    GFED_name=f'/GFED5.1_daily_'
    GFED_NRT_name=f'GFED5NRTspe_CMB_'

if version=='GFED4':
    import h5py
    GFED_path='/work/bb1170/RUN/b383736/data/GFED/GFED4/'
    GFED_name=f'/GFED5.1_daily_'
    


#spath='/work/bb1170/RUN/b383736/data/GFED/GFED5.1_co2_cut_2x2.nc'
spath='/work/bb1170/RUN/b383736/data/GFED/GFED5_co2_larger_2x2_2008_2022.nc'


ds=[]
if version=='GFED5':
    for month in pd.date_range(get_start_date_of_week(start_date)-dt.timedelta(days=10), end_date+dt.timedelta(days=45), freq='M'):
        print(f'{month=}')
        if month.year==2023:
            ds_temp=[]
            for day in pd.date_range(dt.date(month.year,month.month,1),dt.date(month.year,month.month,31)):
                temp=xr.open_dataset(GFED_NRT_path+GFED_NRT_name+f'{month.year}-{month.month:02d}-{day.day:02d}.nc')
                ds_temp.append(temp)
            temp=xr.concat(ds_temp, dim='time', data_vars='all', join='outer')
        else:
            temp=xr.open_dataset(GFED_path+GFED_name+f'{month.year}-{month.month:02d}.nc')
        temp=temp.where((temp.lat>=analy_region_lat[0])&(temp.lat<=analy_region_lat[1])&(temp.lon>=analy_region_lon[0])&(temp.lon<=analy_region_lon[1]),drop=True)
        temp=temp.drop_vars([v for v in temp.data_vars if v not in vars])
        temp['CO2']=(temp.CO2*1e-3/(24*60*60)).assign_attrs(units='kgCO2/s')
        ds.append(temp)
    data=xr.concat(ds, dim='time', data_vars='all', join='outer')

if version=='GFED4':
    year=''
    ds=[]
    for month in pd.date_range(get_start_date_of_week(start_date)-dt.timedelta(days=45), end_date+dt.timedelta(days=10), freq='MS'):
        print(f'{month=}')
        if month.year!=year:
            year=month.year
            temp = h5py.File(f'{GFED_path}/GFED4.1s_{year}_beta.hdf5', 'r')
        ds.append(temp['emissions'][f'{month.month:02d}']['C'][:])
    area = temp['ancill']['grid_cell_area'][:]
    lat = temp['lat'][:]
    lon = temp['lon'][:]
    temp.close()
    C = np.stack(ds)
    data=xr.Dataset({
        'CO2': (['time','lat','lon'],C), 
        'grid_cell_area': (['lat','lon'],area),},
        coords={
            'time': pd.date_range(get_start_date_of_week(start_date)-dt.timedelta(days=45), end_date+dt.timedelta(days=10), freq='MS'), 
            'lat':(['lat', 'lon'], lat),
            'lon':(['lat', 'lon'], lon)})
    data['lat']=data.lat.isel(lon=0)
    data['lon']=data.lon.isel(lat=0)
    data=data.where((data.lat>=analy_region_lat[0])&(data.lat<=analy_region_lat[1])&(data.lon>=analy_region_lon[0])&(data.lon<=analy_region_lon[1]),drop=True)
    days = data.time.dt.days_in_month
    data['CO2']=(data.CO2*1e-3/(days*24*60*60)).assign_attrs(units='kgCO2/m^2s')
    # Convert monthly data to daily by forward-filling values
    data = data.reindex(time=pd.date_range(min(data.time.dt.date.values) , max(data.time.dt.date.values)+dt.timedelta(days=int(days.values[-1])-1) ), method='ffill')
# get weekly sum
# Create the 7-day period bins
period_bins = pd.date_range(start=get_start_date_of_week(start_date), end=(get_start_date_of_week(end_date))+dt.timedelta(days=7), freq="7D")
# Cut the time array into the defined 7-day bins
time_bins = pd.cut(data["time"], bins=period_bins, right=False, labels=pd.to_datetime(period_bins[:-1]))
# Assign the new time bins as coordinates
data=data.assign_coords(time=time_bins)
# weekly sum
# the timestamp corresponds to the first day of the week
data=data.groupby("time").sum('time')
data=data.assign_coords(time=data.time.astype("datetime64[ns]"))
data=data.where((data.time>=pd.to_datetime(start_date-dt.timedelta(days=6)))&(data.time<=pd.to_datetime(end_date)),drop=True)
data['latitude']=data.lat
data['longitude']=data.lon
data=data.swap_dims(lat='latitude', lon='longitude')
data=data.drop_vars(['lat','lon'])
if version=='GFED4':
    # area weighted mean over coarser grid
    data['CO2']=data['CO2']*data.grid_cell_area
    data=data.coarsen(latitude=resolution, longitude=resolution).sum(['latitude','longitude'])
    data['CO2']=(data['CO2']/data.grid_cell_area)
    data.to_netcdf(spath,engine="h5netcdf")
    print(f'saved dataset to {spath}')
elif version=='GFED5':
    data=data.coarsen(latitude=resolution, longitude=resolution).sum(['latitude','longitude'])
    data.to_netcdf(spath)
    print(f'saved dataset to {spath}')

import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt


if __name__ == "__main__":  
    start_date, end_date =dt.date(2010,7,1), dt.date(2010,11,15)
    CAMS_dir='/work/bb1170/RUN/b383736/data/CAMS/satellite/'
    #CAMS_dir='/work/bb1170/RUN/b383736/data/CAMS/IS/'
    region=[2,66,-146,-50]  #lat_min, lat_max, lon_min, lon_max
    k=True
    for day in pd.date_range(start_date,end_date):
        if ((day==day.to_period('M').to_timestamp()) or k):
            if k==True:
                k=False
            else:
                data.close()
            month_str=day.strftime('%Y%m')
            data=xr.open_dataset(f'{CAMS_dir}cams73_v21r2_co2_conc_satellite_inst_{month_str}.nc')
            #data=xr.open_dataset(f'{CAMS_dir}cams73_latest_co2_conc_surface_inst_{month_str}.nc')
            data.where((data.latitude>region[0])& (data.latitude<region[1]),drop=True)
            data.where((data.longitude>region[2])& (data.longitude<region[3]),drop=True)
            #data=data.assign_coords(level=data.level - 1)
            data=data.rename_dims(hlevel='boundaries',time='times',level='levels')
            data = data.rename(time='times',level='levels')
            data['p_boundary']=data.ap+data.bp*data.Psurf
            data['p_diff']=(data.p_boundary.sel(boundaries=slice(0,max(data.boundaries.values)))-data.p_boundary.sel(boundaries=slice(1,max(data.boundaries.values)+1))).rename({'boundaries': 'levels'})
            data['mix']=data.CO2*1e6
            data['xco2']=((data.mix*data.p_diff).sum(dim='levels'))/(data.p_boundary.sel(boundaries=0)-data.p_boundary.sel(boundaries=max(data.boundaries.values)))
            data=data.rename_vars({'Psurf':'pressure'})
        date_str=day.strftime('%Y%m%d')
        temp=data.where((data.times>=day)&(data.times<day+ dt.timedelta(days=1)),drop=True)
        temp.to_netcdf(f"{CAMS_dir}xco2_mean_{date_str}.nc")
        print(f"saving {date_str}" )


import numpy as np
import datetime as dt
import xarray as xr
import pandas as pd
import os

# borders of pressure bins
p_bins=np.array([0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 200, 300, 500, 800, 2000]) # hPa
p_bins=p_bins.dot(100) # hPa -> Pa
num_list=[]
ntot=[]
num_layer=[]
num_tot_layer=[]

# day=20100130
start_date=dt.date(2010,7,1)
# end_date=dt.date(2010,7,3)
end_date=dt.date(2010,8,31)
for d in pd.date_range(start_date,end_date):
# for day in range(20,31):
    print(f'current day: {d}')
    path=f'/work/bb1170/RUN/b382762/data/Flexpart11_output/GOSAT_{d.strftime("%Y_%m")}/Release_{d.strftime("%Y%m%d")}/'
    # path=f'/work/bb1170/RUN/b382762/data/Flexpart11_output/GOSAT_{d.strftime('%Y_%m')}/Release_{d.strftime('%Y%m%d')}/'
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
    
    # particles that reach below 30m sometime
    enhancement_particles=data.where(data.z<30, drop=True).particle.values
    # print(f'day: {day}, {len(data.particle.values)}')
    
    # determine first time
    t_max=pd.to_datetime(np.max(data.time).values)
    # select first time for particles that get below 30m at some point
    data_cut=data.sel(time=t_max,particle=enhancement_particles)
    # data_init=data.sel(time=t_max)  # data for initial particle release # TODO separate different release positions
    # new
    data_init=xr.open_dataset(path+'/part_ic.nc')  # data for initial particle release # TODO separate different release positions
    if data_init.particle.values[0]==0:
        data_init['particle']=data_init.particle+1
    data_init['prs']=p(data_init.height)*100 # p[hPa]->p[Pa]
    
    # number of particles in pressure layer that get below 30m
    num_layer_temp=[len(data_cut.where((data_cut.prs>p_bins[i]) & (data_cut.prs<p_bins[i+1]), drop=True).particle.values) for i in range(0, len(p_bins)-1)]
    # total number of particles released in layer
    num_tot_layer_temp=[len(data_init.where((data_init.prs>p_bins[i]) & (data_init.prs<p_bins[i+1]), drop=True).particle.values) for i in range(0, len(p_bins)-1)]
    
    # number of particles that somewhere reach below 30m
    num=len(data_cut.particle.values)
    num_list.append(num)
    # add total number of particles to list
    ntot.append(len(data_init.particle.values))
    num_layer.append(num_layer_temp)
    num_tot_layer.append(num_tot_layer_temp)
# num_layer

# create dataset
num_ds=xr.Dataset(coords=dict(date=pd.date_range(start_date,end_date), pbin=p_bins[:-1]), 
                  data_vars=dict(num=(['date'],num_list), num_tot=(['date'], ntot), num_layer=(['date', 'pbin'] , num_layer), num_tot_layer=(['date','pbin'], num_tot_layer)))
num_ds['num_layer_rel']=num_ds.num_layer/num_ds.num_tot_layer
num_ds['num_rel']=num_ds.num/num_ds.num_tot
# save dataset
spath='/work/bb1170/RUN/b382762/data/Flexpart11_output/GOSAT_2010_07-08/'
num_ds.to_netcdf(f'{spath}num_enhancement_particles_layers.nc')

# num_ds
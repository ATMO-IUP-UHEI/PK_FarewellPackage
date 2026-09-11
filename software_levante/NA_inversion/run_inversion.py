# run inversions for final flexpart footprints, created 24.03

import xarray as xr
import os
import pandas as pd
import datetime as dt
import numpy as np
import argparse
from pyinverse.loss import Bayesian
from pyinverse.solver import BayesianAnalytical
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
import sys
#sys.path.append("/work/bb1170/RUN/b383736/software/pyinverse/pyinverse/src/pyinverse/")
#from solver import BayesianAnalytical

from utils import get_start_date_of_week, haversine, get_unique_time
import yaml
import argparse
import json
import random
import matplotlib.pyplot as plt
import time
from itertools import product


# get measurements
def get_gosat_measurement_array(start_date, end_date, data_dir,BG="TM5", bg_ds='RemoTeC_2.4.0+IS'):
    '''
    Args:
        start_date: start date of measurements that are used
        end_date: end date of measurements that are used
        data_path: path to the directory containing the TM5-4DVar estimates and backgrounds
        bg_ds: string defining which TM5-4DVar dataset should be used for the background, defaults to RemoTeC_2.4.0+IS
    Returns:
        measurements: array of GOSAT_xco2 - background 
        measurement_covariance: array of GOSAT_xco2_err**2
    
    '''
    # read all measurements
    for date in pd.date_range(start_date, end_date):    
        # check if file exists for that date
        path=f'{data_dir}/{date.strftime("%Y_%m")}/xco2_bg_{date.strftime("%Y%m%d")}.csv'
        if os.path.isfile(path):
            data=pd.read_csv(path)
            if date.date()==start_date:
                gosat_data=pd.DataFrame(data)
            else:
                gosat_data=pd.concat([gosat_data,data])
    gosat_data=gosat_data.reset_index(names='release_num')
    # measurement vector = gosat xco2 - background
    measurements=(gosat_data.xco2-gosat_data[f'{BG}_{bg_ds}_background']).values
    # measurement error
    measurement_covariance=(gosat_data.xco2_err**2).values
    return measurements, measurement_covariance

class measurement_dataset():
    """Measurement dataset for which Flexpart was run

    Args:
        name: name of dataset
        version: version of dataset
        flexpart_path: Flexpart output directory
        flexpart_folder_name: Flexpart subfolder for measurements
        measurement_error: measurement error assumed for all measurements 
        measured: 'co2' (e.g. insitu), 'xco2' (e.g. GOSAT)
        error_inflation: if True, the measurement error gets inflated dependent on the number of other nearby measurements 
        error_inflation_length: distance for measurements to be considered near [km]
        error_inflation_time: time difference for measurements to be considered near [hours]
    """
    
    def __init__(
            self, 
            name: str, 
            version: str, 
            flexpart_path: str,
            flexpart_folder_name: str,
            measurement_error,
            measured: str,
            error_inflation: bool = False,
            error_inflation_length = '',
            error_inflation_time = ''
            ):
        
        self.name, self.version = name, version
        self.flexpart_path = flexpart_path
        self.flexpart_folder_name = flexpart_folder_name
        self.measurement_error = measurement_error
        self.measured = measured
        self.error_inflation = error_inflation
        if error_inflation:
            self.error_inflation_length = error_inflation_length
            self.error_inflation_time = error_inflation_time
    def get_dimension(self,data):
        '''
        Saves size of dataset in self.size
        '''
        self.size=len(data.pointspec.values)
    def get_measurement_error_inflation_array(self,measurement_error, data):
        N=[]
        for i in range(len(data.pointspec.values)):
            measurement=data.isel(pointspec=i)
            temp=data.where(np.abs((data.release_time-measurement.release_time))<=self.error_inflation_time,drop=True)
            dist=haversine(temp.release_lat,temp.release_lon, measurement.release_lat,measurement.release_lon)
            N.append(len(temp.xco2.where(dist<=self.error_inflation_length,drop=True)))
        measurement_error_inflated=np.sqrt(np.array(N))*measurement_error
        return measurement_error_inflated
    
        
def get_is_measurement_array(start_date, end_date, data_dir,BG="TM5", bg_ds='RemoTeC_2.4.0+IS'):
    '''
    Args:
        start_date: start date of measurements that are used
        end_date: end date of measurements that are used
        data_dir: path to the directory containing the Insitu measurements and backgrounds
        bg_ds: string defining which TM5-4DVar dataset should be used for the background, defaults to 'RemoTeC_2.4.0+IS'
    Returns:
        measurements: array of insitu_co2 - background 
    
    '''
    # read all measurements
    data_list=[]
    for date in pd.date_range(start_date, end_date):    
        # check if file exists for that date
        path=f'{data_dir}/{date.strftime("%Y_%m")}/co2_bg_{date.strftime("%Y%m%d")}.csv'
        if os.path.isfile(path):
            data=pd.read_csv(path)
            data_list.append(data)
    is_data=pd.concat(data_list)
    is_data=is_data.reset_index(names='release_num')
    # measurement vector = gosat xco2 - background
    measurements=(is_data['co2_val[ppm]']-is_data[f'{BG}_{bg_ds}_background']).values
    return measurements
# get prior and prior covariance
def get_weekly_priors_from_flux(prior_flux_path, start_date, end_date):
    ''' gets weekly TM5-4DVar prior and flat prior (area weighted mean, for month transition, weighted mean for number of days in each month)
        from prior_flux_path dataset for desired time period
    Args:
        prior_flux_path: path to '...{res}x{res}_cut.nc' flux dataset, with res spatial resolution of statevector
        start_date: start date of measurements that are used
        end_date: end date of measurements that are used
    Returns:
        flat_prior: array of length of spatial res*number of weeks, with mean flux for each week in selected time period
        TM5_prior: array of length of spatial res*number of weeks, with TM5-4DVar flux for each week
        prior_flux_sel: weekly prior
    '''
    # footprints start 10 days before startdate
    f_start=start_date-dt.timedelta(days=10)
    # get first day of the week containing the start of footprints
    f_start=get_start_date_of_week(f_start)
    # Create the 7-day period bins
    period_bins = pd.date_range(start=f_start, end=end_date+dt.timedelta(days=7), freq="7D")

    # read prior
    prior_flux=xr.open_dataset(prior_flux_path)
    # select months
    prior_flux_sel=prior_flux.sel(time=slice(f"{f_start.year}-{(f_start.month)}",f"{end_date.year}-{end_date.month}"))

    # Convert monthly data to daily by forward-filling values
    prior_flux_sel = prior_flux_sel.reindex(time=pd.date_range(f_start,end_date), method='ffill')   #  .sel(latitude=9.5, longitude=-136.5).total_flux.values
    # Cut the time array into the defined 7-day bins
    time_bins = pd.cut(prior_flux_sel.time, bins=period_bins, right=False, labels=period_bins[:-1],include_lowest=True)
    # Assign the new time bins as coordinates
    prior_flux_sel=prior_flux_sel.assign_coords(time=time_bins)
    # weekly mean, only different value for weeks across two months
    # the timestamp corresponds to the first day of the week
    prior_flux_sel=prior_flux_sel.groupby("time").mean('time')

    # area weighted mean for flat prior
    # flux in kg/(m^2 s)
    mean_prior=((prior_flux_sel.total_flux*prior_flux_sel.grid_cell_area).sum(dim=['latitude','longitude'])/prior_flux_sel.grid_cell_area.sum(dim=['latitude','longitude']))
    flat_prior=((np.ones([(prior_flux_sel.sizes['latitude']*prior_flux_sel.sizes['longitude']),prior_flux_sel.sizes['time']])*mean_prior.values).T).flatten()
    
    # use TM5-4DVar prior flux
    TM5_prior=prior_flux_sel.total_flux.values.flatten()
    
    return flat_prior, TM5_prior, prior_flux_sel
def get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=3, epsilon=0.84, prior_min=0.005,corr='both'):
    ''' 
    Args:
        weekly_prior: prior fluxes, with dimesions time, latitude,longitude
        L: covariance lenght parameter, defaults to 500 km
        T: covariance time parameter, defaults to 3 months
        epsilon: fraction of prior flux used for prior_std, defaults to 0.84
        prior_min: minimum value used for prior_std flux, defaults to 0.005 kg CO2/sec/grid box
    Retruns:
        prior_cov: covariance matrix
    '''
    from scipy import sparse
    if corr=='both' or corr=='temporal':
        T=(T*30.437)    # months -> days
        local_time = time.localtime()
        # get time differces
        t_vals=np.array(weekly_prior.time.values, dtype='datetime64[D]')
        # C_T(t1,t2)= exp(-|t1-t2|/T)
        #C_T=np.exp(-np.abs(t_vals[:, None] - t_vals[None, :])/ np.timedelta64(1, 'D')/T)
        dt = np.abs(t_vals[:, None] - t_vals[None, :]) / np.timedelta64(1, 'D')
        rows,cols = np.where( dt <= 3*T)
        vals = np.exp(-dt[rows,cols]/T)
        C_T = sparse.csr_matrix((vals, (rows, cols)), shape=dt.shape)
        print(C_T.shape)
        del dt
        del vals
        del t_vals

    if corr=='both' or corr=='spatial':
        # get spatial distances
        lat_vals = weekly_prior.latitude.values
        lon_vals = weekly_prior.longitude.values
        
        # Compute geodesic distances vectorized
        #latlon_pairs = np.column_stack([lat_vals, lon_vals])
        latlon_pairs=np.array(list((product(lat_vals,lon_vals))))
        D=squareform(pdist(latlon_pairs, metric=lambda u, v: haversine(u[0], u[1], v[0], v[1])))
        rows,cols = np.where(D<=3*L)
        vals = np.exp(-D[rows,cols]/L)
        C_r = sparse.csr_matrix((vals, (rows, cols)), shape=D.shape)
        del latlon_pairs
        del D
        del rows, cols
        del vals
        #C_r[mask] = np.exp(-D[mask] / L)
        #C_r = np.exp(-squareform(pdist(latlon_pairs,  metric=lambda u, v: dist(u,v))/L))
        
    if corr=='temporal':
        C_r=np.eye(len(weekly_prior.latitude.values)*len(weekly_prior.longitude.values))
    if corr=='spatial':
        C_T=np.eye(len(weekly_prior.time.values))


    # get prior uncertainties
    weekly_prior['prior_min'] =(prior_min/weekly_prior.grid_cell_area).assign_attrs(units='kg_CO2/(m^2 s)') # 1/grid_box = 1/a 1/m^2
    # returns max val of weekly_prior.total_flux*epsilon and weekly_prior.prior_min
    weekly_prior['prior_std'] =  xr.where((np.abs(weekly_prior.total_flux*epsilon)>np.abs(weekly_prior.prior_min)), np.abs(weekly_prior.total_flux*epsilon), np.abs(weekly_prior.prior_min))
    # prior covariance
    # cov (x_r1,t1, x_r2,t2)=sig_r1,t1 * sig_r1,t2 * C_r(r1,r2) * C_T(t1,t2)
    # prior_var = sig_r1,t1 * sig_r1,t2
    weekly_prior=weekly_prior.stack(grid_box=("time","latitude", "longitude")).squeeze()
    prior_var=((np.array(weekly_prior['prior_std'].values)*(np.array(weekly_prior['prior_std'].values).T)))
    prior_var=np.eye(len(prior_var))*prior_var
    #prior_cov=np.matmul(prior_var,np.multiply(C_r, C_T)) 
    # matmul: matrix multiplication
    # multiply elementwise multiplication 
    # total correlation C(x_r1,t1, x_r2,t2)= C_r(r1,t1,r2,t2) * C_T(r1,t1,r2,t2) (elementwise multiplication)
    # covariance: cov (x_r1,t1, x_r2,t2)=sig_r1,t1  * C(r1,t1,r2,t2)* sig_r2,t2 = Sig * C * Sig

    sig=sparse.diags(weekly_prior.prior_std.values)
    

    C = sparse.kron(C_T, C_r, format="csr")
    prior_cov=sig @ C @sig
    del sig
    local_time = time.localtime()
    print("Aktuelle Uhrzeit:", time.strftime("%d-%H:%M:%S", local_time))
    if False: # save correlation matrices?
        C_r = xr.DataArray(C_r.todense(), dims=['x', 'y'])
        C_r.to_netcdf(f'/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/inversions_2M/2x2/with_fixed_correlation/corr_cr_{corr}_corr_2x2.nc')
        C_T = xr.DataArray(C_T.todense(), dims=['x', 'y'])
        C_T.to_netcdf(f'/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/inversions_2M/2x2/with_fixed_correlation/corr_ct_{corr}_corr_2x2.nc')
    del C_r
    del C_T
    return prior_cov

def get_cov_from_weekly_prior(weekly_prior, L=500,T=3, epsilon=0.84, prior_min=0.005,corr='both'):
    ''' 
    Args:
        weekly_prior: prior fluxes, with dimesion already stacked (.stack(grid_box=("time","latitude", "longitude")).squeeze())
        L: covariance lenght parameter, defaults to 500 km
        T: covariance time parameter, defaults to 3 months
        epsilon: fraction of prior flux used for prior_std, defaults to 0.84
        prior_min: minimum value used for prior_std flux, defaults to 0.005 kg CO2/sec/grid box
    Retruns:
        prior_cov: covariance matrix
    '''
    local_time = time.localtime()
    print("Aktuelle Uhrzeit:", time.strftime("%H:%M:%S", local_time))
    if corr=='with_fixed':
        T=(T*30.437)    # months -> days
        local_time = time.localtime()
        print("Aktuelle Uhrzeit:", time.strftime("%H:%M:%S", local_time))
        # get time differces
        t_vals=np.array(weekly_prior.time.values, dtype='datetime64[D]')
        # C_T(t1,t2)= exp(-|t1-t2|/T)
        #C_T=np.exp(-np.abs(t_vals[:, None] - t_vals[None, :])/ np.timedelta64(1, 'D')/T)
        dt = np.abs(t_vals[:, None] - t_vals[None, :]) / np.timedelta64(1, 'D')
        rows,cols = np.where( dt <= 2*T)
        vals = np.exp(-dt[rows,cols]/T)
        C_T = sparse.csr_matrix((vals, (rows, cols)), shape=dt.shape)
        print(C_T.shape)
        del dt
        del vals
        del t_vals
        # get spatial distances
        lat_vals = weekly_prior.latitude.values
        lon_vals = weekly_prior.longitude.values
        dlat_km=haversine(np.mean(lat_vals),np.mean(lon_vals),np.mean(lat_vals)+res,np.mean(lon_vals))/res
        dlon_km=haversine(np.mean(lat_vals),np.mean(lon_vals),np.mean(lat_vals),np.mean(lon_vals)+res)/res
        # Compute geodesic distances vectorized
        latlon_pairs = np.column_stack([lat_vals, lon_vals])
        # C_r(r1,r2)= exp(-|r1-r2|/L)
        def dist(u, v ):
            return np.sqrt((u[0]-v[0])**2*dlat_km**2 + (u[1]-v[1])**2*dlon_km**2)
        D=squareform(pdist(latlon_pairs,  metric=lambda u, v: dist(u,v)))
        rows,cols = np.where(D<=2*L)
        vals = np.exp(-D[rows,cols]/L)
        C_r = sparse.csr_matrix((vals, (rows, cols)), shape=D.shape)
        del latlon_pairs
        del D
        del rows, cols
        del vals
        #C_r[mask] = np.exp(-D[mask] / L)
        #C_r = np.exp(-squareform(pdist(latlon_pairs,  metric=lambda u, v: dist(u,v))/L))
        

    if corr=='both' or corr=='temporal':
        T=(T*30.437)    # months -> days
        # get time differces
        t_vals=np.array(weekly_prior.time.values, dtype='datetime64[D]')
        # C_T(t1,t2)= exp(-|t1-t2|/T)
        C_T=np.exp(-np.abs(t_vals[:, None] - t_vals[None, :])/ np.timedelta64(1, 'D')/T)
    if corr=='both' or corr=='spatial':
        # get spatial distances
        lat_vals = weekly_prior.latitude.values
        lon_vals = weekly_prior.longitude.values
        # Compute geodesic distances vectorized
        latlon_pairs = np.column_stack([lat_vals, lon_vals])
        # C_r(r1,r2)= exp(-|r1-r2|/L)
        C_r = np.exp(-squareform(pdist(latlon_pairs, metric=lambda u, v: haversine(u[0], u[1], v[0], v[1])))/L)
    # get prior uncertainties
    weekly_prior['prior_min'] =(prior_min/weekly_prior.grid_cell_area).assign_attrs(units='kg_CO2/(m^2 s)') # 1/grid_box = 1/a 1/m^2
    # returns max val of weekly_prior.total_flux*epsilon and weekly_prior.prior_min
    weekly_prior['prior_std'] =  xr.where((np.abs(weekly_prior.total_flux*epsilon)>np.abs(weekly_prior.prior_min)), np.abs(weekly_prior.total_flux*epsilon), np.abs(weekly_prior.prior_min))
    # prior covariance
    # cov (x_r1,t1, x_r2,t2)=sig_r1,t1 * sig_r1,t2 * C_r(r1,r2) * C_T(t1,t2)
    # prior_var = sig_r1,t1 * sig_r1,t2
    prior_var=((np.array(weekly_prior['prior_std'].values)*(np.array(weekly_prior['prior_std'].values).T)))
    prior_var=np.eye(len(prior_var))*prior_var
    #prior_cov=np.matmul(prior_var,np.multiply(C_r, C_T)) 
    # matmul: matrix multiplication
    # multiply elementwise multiplication 
    # total correlation C(x_r1,t1, x_r2,t2)= C_r(r1,t1,r2,t2) * C_T(r1,t1,r2,t2) (elementwise multiplication)
    # covariance: cov (x_r1,t1, x_r2,t2)=sig_r1,t1  * C(r1,t1,r2,t2)* sig_r2,t2 = Sig * C * Sig
    if corr=='both':
        prior_cov=np.matmul(np.diag(weekly_prior.prior_std.values),np.matmul(np.multiply(C_r,C_T), np.diag(weekly_prior.prior_std.values)))
    if corr=='temporal':
        prior_cov=np.matmul(np.diag(weekly_prior.prior_std.values),np.matmul(C_T, np.diag(weekly_prior.prior_std.values)))
    if corr=='spatial':
        prior_cov=np.matmul(np.diag(weekly_prior.prior_std.values),np.matmul(C_r, np.diag(weekly_prior.prior_std.values)))
    if corr=='with_fixed':
        #prior_cov=np.matmul(np.diag(weekly_prior.prior_std.values),np.matmul(np.multiply(C_r,C_T), np.diag(weekly_prior.prior_std.values)))
        sig=sparse.diags(weekly_prior.prior_std.values)
        #prior_cov= C_r.multiply(C_T)
        #prior_cov=prior_cov.multiply(weekly_prior.prior_std.values[:,None])
        #prior_cov=prior_cov.multiply(weekly_prior.prior_std.values[None,:])
        prior_cov=sig @ (C_r.multiply(C_T))@sig
        del sig
        #prior_cov=weekly_prior.prior_std.values[:,None]* C_r.multiply(C_T)* weekly_prior.prior_std.values[None,:]
    local_time = time.localtime()
    print("Aktuelle Uhrzeit:", time.strftime("%H:%M:%S", local_time))
    return prior_cov.toarray()

def get_prior_var_no_correlation_from_weekly_prior(weekly_prior, epsilon=0.84, prior_min=0.005):
    ''' Get prior variance from weekly prior fluxes, no correlation 
    Args:
        weekly_prior: prior fluxes, with dimesion already stacked (.stack(grid_box=("time","latitude", "longitude")).squeeze())
        epsilon: fraction of prior flux used for prior_std, defaults to 0.84
        prior_min: minimum value used for prior_std flux, defaults to 0.005 kg CO2/sec/grid box
    Retruns:
        prior_var: covariance matrix
    '''
    # get prior uncertainties
    weekly_prior['prior_min'] =(prior_min/weekly_prior.grid_cell_area).assign_attrs(units='kg_CO2/(m^2 s)') # 1/grid_box = 1/a 1/m^2
    # returns max val of weekly_prior.total_flux*epsilon and weekly_prior.prior_min
    weekly_prior['prior_std'] =  xr.where((np.abs(weekly_prior.total_flux*epsilon)>np.abs(weekly_prior.prior_min)), np.abs(weekly_prior.total_flux*epsilon), np.abs(weekly_prior.prior_min))
    prior_var=sparse.diags((weekly_prior['prior_std'].values.flatten()**2), format='csr')
    return prior_var

def run_inv_TM5_prior_flat(flat_prior, TM5_4DVar_prior, prior_covariance, measurements,measurement_covariance, footprint,spath,SAVE_AK=False, footprint_col_name='spec001_mr'):
    ''' function that runs two inversions (two different priors), saves output as dataset
    Args:
        flat_prior: flat_prior as array
        TM5_4DVar_prior: TM5_4DVar_prior as array
        prior_covariance: prior_covariance as array, used for both priors
        measurements: measurements as array
        measurement_covariance: measurement_covariance as array
        footprint: footprint as dataarray with dimension 'grid_box' = Multiindex of time, latitude and longitude
        footprint_col_name: footprint column name, defaults to 'spec001_mr'
        spath: path where output dataset is saved
        SAVE_AK: set to True is temporal mean of AK and posterior cov should be saved
    Returns:
        nothing, saves dataset
    '''
    
    # inversion with flat prior
    flat_loss = Bayesian(
        x_prior=flat_prior,
        cov_prior=prior_covariance,
        y=measurements,
        cov_y=measurement_covariance,
        K=footprint[footprint_col_name].values,
    )
    flat_solver = BayesianAnalytical(flat_loss)

    # inversion with TM5-4DVar prior
    TM5_loss = Bayesian(
        x_prior=TM5_4DVar_prior,
        cov_prior=prior_covariance,
        y=measurements,
        cov_y=measurement_covariance,
        K=footprint[footprint_col_name].values,
    )
    TM5_solver = BayesianAnalytical(TM5_loss)
    # add data to ds
    ds=xr.Dataset(data_vars=dict(
            flat_prior_flux=(["grid_box"], flat_prior,{"units": "kgCO2/(m^2 s)"}),
            TM5_prior_flux=(["grid_box"], TM5_4DVar_prior,{"units": "kgCO2/(m^2 s)"}),
            prior_uncertainty=(["grid_box"], np.sqrt(prior_covariance.diagonal()),{"units": "kgCO2/(m^2 s)"}),
            # flat prior
            flat_posterior_flux=(["grid_box"], flat_solver.x_posterior,{"units": "kgCO2/(m^2 s)"}),
            flat_posterior_std=(["grid_box"], np.sqrt(np.diag(flat_solver.cov_posterior)),{"units": "kgCO2/(m^2 s)"}),
            flat_averaging_kernel_diag=(["grid_box"], np.diag(flat_solver.averaging_kernel),{"units": "kgCO2/(m^2 s)"}),
            # TM5-4DVar prior
            TM5_posterior_flux=(["grid_box"], TM5_solver.x_posterior,{"units": "kgCO2/(m^2 s)"}),
            TM5_posterior_std=(["grid_box"], np.sqrt(np.diag(TM5_solver.cov_posterior)),{"units": "kgCO2/(m^2 s)"}),
            TM5_averaging_kernel_diag=(["grid_box"], np.diag(TM5_solver.averaging_kernel),{"units": "kgCO2/(m^2 s)"}),
            # measurements
            meas=(['meas_num'], measurements, {"units": "ppm"}),
            meas_cov=(['meas_num'], measurement_covariance, {"units": "ppm^2"}),
        ),
        coords=dict(
            grid_box=footprint.grid_box,
            meas_num=np.arange(0,measurements.size,step=1),
        ))
    ds=ds.unstack(dim='grid_box')
    # ds[['flat_prior_flux','flat_posterior_flux','TM5_prior_flux','TM5_posterior_flux','prior_uncertainty']] = ds[['flat_prior_flux','flat_posterior_flux','TM5_prior_flux','TM5_posterior_flux','prior_uncertainty']].assign_attrs(units='')
    ds.to_netcdf(spath)
    print(f'saved dataset to: {spath}')
    
    if SAVE_AK:
        # get temporal mean of ak and posterior cov
        ak_ds=xr.Dataset(data_vars=dict(
            # necessary for flat prior inversion?
            # flat_averaging_kernel=(["x","y"], flat_solver.averaging_kernel),
            TM5_posterior_cov=(["x","y"], TM5_solver.cov_posterior),
            TM5_averaging_kernel=(["x","y"], TM5_solver.averaging_kernel)
            ))
        # First, copy the coordinates from footprints to ak_ds
        # ak_ds = ak_ds.assign_coords(
        #     x=footprints.grid_box,
        #     y=footprints.grid_box)

        # # Now, 'grid_box' is available; we can unstack
        # ak_ds = ak_ds.swap_dims({"x": "grid_box", "y": "grid_box"})  # if x and y both are 1D parts of grid_box

        # # Now unstack 'grid_box' into 'time', 'latitude', 'longitude'
        # ak_ds = ak_ds.unstack("grid_box")

        # # Now you can average over 'time'
        # ak_ds_mean = ak_ds.mean(dim="time")
        spath_mean=spath.replace('.nc', '_ak_post_cov.nc')
        print(f'saving to {spath_mean}')
        ak_ds.to_netcdf(spath_mean)
        print('successfull')
    return

def run_inv(TM5_4DVar_prior, prior_covariance, measurements,measurement_covariance, footprint,spath,SAVE_AK=False, footprint_col_name='spec001_mr'):
    ''' function that runs inversion, saves output as dataset
    Args:
        TM5_4DVar_prior: TM5_4DVar_prior as array
        prior_covariance: prior_covariance as array, used for both priors
        measurements: measurements as array
        measurement_covariance: measurement_covariance as array
        footprint: footprint as dataarray with dimension 'grid_box' = Multiindex of time, latitude and longitude
        footprint_col_name: footprint column name, defaults to 'spec001_mr'
        spath: path where output dataset is saved
        SAVE_AK: set to True is temporal mean of AK and posterior cov should be saved
    Returns:
        nothing, saves dataset
    '''
    print('starting inversion')
    # inversion with TM5-4DVar prior
    TM5_loss = Bayesian(
        x_prior=TM5_4DVar_prior,
        cov_prior=prior_covariance,
        y=measurements,
        cov_y=measurement_covariance,
        K=footprint[footprint_col_name].values,
    )
    print('calculated Loss')
    TM5_solver = BayesianAnalytical(TM5_loss)
    print(prior_covariance.diagonal().shape)
    local_time = time.localtime()
    print("Aktuelle Uhrzeit:", time.strftime("%d-%H:%M:%S", local_time))
    # add data to ds
    ds=xr.Dataset(data_vars=dict(
            TM5_prior_flux=(["grid_box"], TM5_4DVar_prior,{"units": "kgCO2/(m^2 s)"}),
            prior_uncertainty=(["grid_box"], np.sqrt(prior_covariance.diagonal()),{"units": "kgCO2/(m^2 s)"}),
            # TM5-4DVar prior
            TM5_posterior_flux=(["grid_box"], TM5_solver.x_posterior,{"units": "kgCO2/(m^2 s)"}),
            TM5_posterior_std=(["grid_box"], TM5_solver.std_posterior,{"units": "kgCO2/(m^2 s)"}),
            TM5_averaging_kernel_diag=(["grid_box"], TM5_solver.averaging_kernel.diagonal(),{"units": "kgCO2/(m^2 s)"}),
            # measurements
            meas=(['meas_num'], measurements, {"units": "ppm"}),
            meas_cov=(['meas_num'], measurement_covariance, {"units": "ppm^2"}),
        ),
        coords=dict(
            grid_box=footprint.grid_box,
            meas_num=np.arange(0,measurements.size,step=1),
        ))
    local_time = time.localtime()
    print("Aktuelle Uhrzeit:", time.strftime("%d-%H:%M:%S", local_time))
    ds=ds.unstack(dim='grid_box')
    ds.to_netcdf(spath)
    print(f'saved dataset to: {spath}')
    if SAVE_AK:
        # get temporal mean of ak and posterior cov
        ak_ds=xr.Dataset(data_vars=dict(
            # necessary for flat prior inversion?
            # flat_averaging_kernel=(["x","y"], flat_solver.averaging_kernel),
            TM5_posterior_cov=(["x","y"], TM5_solver.cov_posterior),
            TM5_averaging_kernel=(["x","y"], TM5_solver.averaging_kernel)
            ))
        # First, copy the coordinates from footprints to ak_ds
        # ak_ds = ak_ds.assign_coords(
        #     x=footprints.grid_box,
        #     y=footprints.grid_box)

        # # Now, 'grid_box' is available; we can unstack
        # ak_ds = ak_ds.swap_dims({"x": "grid_box", "y": "grid_box"})  # if x and y both are 1D parts of grid_box

        # # Now unstack 'grid_box' into 'time', 'latitude', 'longitude'
        # ak_ds = ak_ds.unstack("grid_box")

        # # Now you can average over 'time'
        # ak_ds_mean = ak_ds.mean(dim="time")
        spath_mean=spath.replace('.nc', '_ak_post_cov.nc')
        print(f'saving to {spath_mean}')
        ak_ds.to_netcdf(spath_mean)
        print('successfull')
    return

def plot_prior(prior,flattened=False,prior_not_flattened=None,name='prior_flux.png',spath='/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/figures/'):
    if flattened:
        dims=['time','latitude','longitude']
        coords = {
            'time': prior_not_flattened.time.values,  
            'latitude': prior_not_flattened.latitude.values,
            'longitude': prior_not_flattened.longitude.values
        }
        prior=xr.DataArray(prior.reshape((len(prior_not_flattened.time.values),len(prior_not_flattened.latitude.values),len(prior_not_flattened.longitude.values))), dims=dims,coords=coords)
    analy_region_lat=[20,48]
    analy_region_lon=[-126,-70]

    res=2
    n_regions=[2,3] #lat, lon
    lat_dist=int((analy_region_lat[1]-analy_region_lat[0])/n_regions[0]/res)*res
    lon_dist=round((analy_region_lon[1]-analy_region_lon[0])/n_regions[1]/res)*res
    #define region edges
    regions_lat=[]
    for i in range(n_regions[0]+1): #lat
        if i == n_regions[0]:
            regions_lat.append(analy_region_lat[1])
        else:
            regions_lat.append(analy_region_lat[0]+i*lat_dist)
    regions_lon=[]
    for i in range(n_regions[1]+1): #lat
        if i == n_regions[1]:
            regions_lon.append(analy_region_lon[1])
        else:
            regions_lon.append(analy_region_lon[0]+i*lon_dist)
    fig,axs=plt.subplots(2,3,figsize=(24,12))
    for i in range(n_regions[0]): #lat
        for j in range(n_regions[1]): #lon
            k=0

            temp3=prior.where(((prior.latitude>= regions_lat[i]) &(prior.latitude<regions_lat[i+1])),drop=True)
            temp3=temp3.where(((temp3.longitude>= regions_lon[j]) &(temp3.longitude<regions_lon[j+1])),drop=True)
            temp3=temp3.mean(['latitude','longitude'])
            temp3['total_flux']=temp3*1e6
            #temp3['prior_uncertainty']=temp3.prior_uncertainty*1e6
            axs[i,j].plot(temp3.time, temp3, c='orange', marker='o', markersize=2, label='prior TM5-4DVar flux')


            #plt.fill_between(temp3.time, temp3.TM5_posterior_flux-temp2.TM5_posterior_std, temp2.TM5_posterior_flux+temp2.TM5_posterior_std, color='red', alpha=0.2) 
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{spath}{name}',bbox_inches='tight')
    plt.close()
# utils
def load_config(config_path):
    """Load configuration with from config YAML file."""
    with open(config_path, 'r') as file:
        raw_config = yaml.safe_load(file)
    # Parse dates
    raw_config['start_date'] = dt.datetime.strptime(raw_config['start_date'], "%Y-%m-%d").date()
    raw_config['end_date'] = dt.datetime.strptime(raw_config['end_date'], "%Y-%m-%d").date()
    return raw_config

# pass config file path from slurm skript
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, help='Path to YAML config file')
    return parser.parse_args()


if __name__ == "__main__":
    # get config file path
    args = parse_args()
    config_path = args.config
    
    # read config from file
    CONFIG = load_config(config_path)
    for key, value in CONFIG.items():
        globals()[key] = value
    print(inversion_subdirectory)

    gosat=measurement_dataset(
        name = 'gosat', 
        version = '2.4.1', 
        flexpart_path = output_dir+'/',
        flexpart_folder_name = 'RemoTeCv240',
        measurement_error=gosat_meas_err_list,
        measured = 'xco2',
        error_inflation = True,
        error_inflation_length=500,
        error_inflation_time=np.timedelta64(dt.timedelta(hours=1))
        )
    
    insitu=measurement_dataset(
        name='insitu',
        version='GLOBALVIEWplus_v10.1',
        flexpart_path = output_dir+'/',
        flexpart_folder_name = 'insitu',
        measurement_error=insitu_meas_err_list,
        measured = 'co2',
        error_inflation = False,)
    
    datasets= [gosat,insitu]

    for res in res_list:             # ,2, 4
        # read data
        measurements=[]
        data=[]
        for dataset in datasets:
            dataset_path=f'{dataset.flexpart_path}/{dataset.flexpart_folder_name}/{footprints_subdirectory}/high_res_scaled_footprints_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{res}x{res}_weekly.nc'
            data_temp= xr.open_dataset(dataset_path)
            # cut to smaller area
            if SMALER_MEAS_AREA:
                data_temp=data_temp.where((data_temp.release_lat>interp_region[0])&(data_temp.release_lat<interp_region[1]) & 
                                            (data_temp.release_lon>interp_region[2])&(data_temp.release_lon<interp_region[3]), drop=True)
            if FILTER_GOSAT_MEAS: # filter, only use GOSAT meas, if more than specific value per week
                if dataset.name=='gosat':
                    data_temp['release_time'] = xr.apply_ufunc(get_unique_time,data_temp['release_time'],input_core_dims=[['time']],vectorize=True,dask='parallelized')
                    data_temp['release_day'] = xr.apply_ufunc(get_unique_time,data_temp['release_day'],input_core_dims=[['time']],vectorize=True,dask='parallelized')
                    data_temp['pointspec']=np.arange(0,data_temp.pointspec.size)    # assign integer to each measurement
                    # get release week and release box
                    data_temp['release_week']=(('pointspec'), [get_start_date_of_week(pd.to_datetime(data_temp['release_day'][i].item())) for i in range(data_temp.pointspec.size)])
                    data_temp['release_box_lat']=(('pointspec', ), [data_temp.latitude[np.abs(data_temp.release_lat.values[i]-data_temp.latitude).argmin()].item() for i in range(data_temp.pointspec.size)])
                    data_temp['release_box_lon']=(('pointspec', ), [data_temp.longitude[np.abs(data_temp.release_lon.values[i]-data_temp.longitude).argmin()].item() for i in range(data_temp.pointspec.size)])
                    # create dataframe
                    df = xr.Dataset({
                        'release_week': data_temp.release_week,
                        'release_box_lat': data_temp.release_box_lat,
                        'release_box_lon': data_temp.release_box_lon
                    }).to_dataframe().reset_index()
                    # Count occurrences of each (release_week, release_box) group
                    counts = df.groupby(['release_week', 'release_box_lat', 'release_box_lon']).size().reset_index(name='count')
                    # Filter to keep only duplicates (count > min_num_gosat_meas, defined above)
                    duplicates = counts[counts['count'] > min_num_gosat_meas][['release_week', 'release_box_lat', 'release_box_lon']]
                    # Merge to find matching pointspecs
                    merged = df.merge(duplicates, on=['release_week', 'release_box_lat', 'release_box_lon'])
                    # Extract the pointspec indices to keep
                    pointspec_to_keep = merged['pointspec'].unique()
                    # filter dataset
                    measurements.extend( data_temp.sel(pointspec=pointspec_to_keep)  )          
                    # same for with and without correlation
            # get measurements
            elif VERIFICATION_SAMPLE:
                if os.path.exists(f'{output_dir}{inversion_subdirectory}/verification.json'):
                    with open(f'{output_dir}{inversion_subdirectory}/verification.json', 'r', encoding='utf-8') as file:
                        indizes=json.load(file)
                else:
                    verification_sample_gosat=random.sample(list(range(0,len(gosat_data.pointspec))),round(len(gosat_data.pointspec)*verification_gosat))
                    dataset_sample_gosat=list(range(0,len(gosat_data.pointspec)))
                    verification_sample_insitu=random.sample(list(range(0,len(is_data.pointspec))),round(len(is_data.pointspec)*verification_insitu))
                    dataset_sample_insitu=list(range(0,len(is_data.pointspec)))
                    for s in verification_sample_gosat:
                        dataset_sample_gosat.remove(s)
                    for s in verification_sample_insitu:
                        dataset_sample_insitu.remove(s)
                    indizes={ "verification_indizes_gosat": verification_sample_gosat , "dataset_indizes_gosat":dataset_sample_gosat,
                                "verification_indizes_insitu":verification_sample_insitu, "dataset_indizes_insitu":dataset_sample_insitu}
                    if not os.path.exists(f'{output_dir}/{inversion_subdirectory}/'):
                        os.makedirs(f'{output_dir}/{inversion_subdirectory}/')
                    with open(f'{output_dir}/{inversion_subdirectory}/verification.json','w') as f:
                        json.dump(indizes,f)
                        print(f'saved indizes of verification sample in {output_dir}/{inversion_subdirectory}/verification.json')
                data_temp = data_temp.isel(pointspec=indizes[f'dataset_indizes_{dataset.name}'])
                measurements.extend((data_temp.dataset.measured-data_temp[f'{BG}_{bg_ds}_background']).values)

            else:   # without measurements
                measurements.extend((data_temp[dataset.measured]-data_temp[f'{BG}_{bg_ds}_background']).values)
                dataset.get_dimension(data_temp)

            data.append(data_temp)
        # combine measurement arrays
        measurements=np.array(measurements)
        
        # read prior flux 
        prior_flux_path=f'{prior_flux_dir}/flux_{res}x{res}_prior_cut.nc'
        flat_prior, TM5_prior, weekly_prior=get_weekly_priors_from_flux(prior_flux_path, start_date, end_date)
        # get prior covariance from prior flux
        weekly_prior_flattened=weekly_prior.stack(grid_box=("time","latitude", "longitude")).squeeze()

        # run with / without covariance
        for corr_str in corr_list:     #, 'no'
            print(f'{corr_str} covariance')
            if prior_cov_dir != "":
                cov_path=f'{prior_cov_dir}/'
                if os.path.exists(f'{cov_path}/cov_{corr_str}_corr_{res}x{res}.nc'):
                    print(f'read in covariance matrix from {cov_path}')
                else:
                    cov_path=f'{output_dir}/{inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/'
                    print(f'prior covariance matrix does not exist in directory {cov_path}, checking in inversion directory')
            else:
                cov_path=f'{output_dir}/{inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/'
            if not os.path.isdir(cov_path):
                print(f"make dir: {cov_path}")
                os.makedirs(f"{cov_path}")
            cov_path+=f'cov_{corr_str}_corr_{res}x{res}.nc'
            if os.path.isfile(cov_path):
                prior_cov=xr.open_dataarray(cov_path).values
                print(f'read cov matrix from {cov_path}')
            else:
                if corr_str=='with':
                    print('getting cov matrix')
                    prior_cov = get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=3, epsilon=0.84, prior_min=0.005,corr='both')
                elif corr_str=='with_1M':
                    print('getting cov matrix')
                    prior_cov = get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=1, epsilon=0.84, prior_min=0.005,corr='both')
                elif corr_str=='with_e04':
                    print('getting cov matrix')
                    prior_cov = get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=3, epsilon=0.4, prior_min=0.005)
                elif corr_str=='no':
                    print('getting prior variance')
                    prior_cov = get_prior_var_no_correlation_from_weekly_prior(weekly_prior)
                    import sparse
                    prior_cov=(sparse.COO.from_scipy_sparse(prior_cov)).todense()
                elif corr_str=='temporal' or corr_str=='spatial':
                    prior_cov = get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=3, epsilon=0.84, prior_min=0.005,corr=corr_str)
                elif corr_str=='temporal_1M':
                    prior_cov = get_cov_from_not_flattened_weekly_prior(weekly_prior, L=500,T=1, epsilon=0.84, prior_min=0.005,corr='temporal')
                else: 
                    print(f'correlation string {corr_str} is not defined')
                
                if False:
                    # save covariance matrix
                    if not os.path.isdir(f'{output_dir}/{inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/'):
                        os.makedirs(f'{output_dir}/{inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/')
                    # create datarray
                    prior_cov_da = xr.DataArray(prior_cov, dims=['x', 'y'])
                    prior_cov_da.to_netcdf(cov_path)
                    print(f'saved covariance matrix to {cov_path}')
                    del prior_cov_da
            print('have prior covariance')
            # run for different footprint scalings
            # only use spec001_mr_scaled_beta_prime with offset
            for f_col in f_list:        
                print(f'current footprint: {f_col}')
                footprints=[]
                measurement_errors=[]
                for i, dataset in enumerate(datasets):
                    # get footprints
                    footprints.append(data[i][[f_col]])
                    measurement_errors.append(dataset.measurement_error)
                # combine footprints, pointspec dim=first all gosat, then all insitu
                footprints=xr.concat(footprints, dim='pointspec')
                footprints=footprints.stack(grid_box=("time","latitude", "longitude")).squeeze()

                measurement_errors=list(product(*measurement_errors))
            
                for measurement_error in measurement_errors:
                    # create directory in which footprints will be saved
                    sdir=f'{output_dir}/{inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_col}/'
                    for i, dataset in enumerate(datasets):
                        print(f'meas_err_val {dataset.name}: {measurement_error[i]}')
                        sdir+=f'{measurement_error[i]}ppm_{dataset.name}_meas_err/'
                        if FILTER_GOSAT_MEAS and dataset.name=='gosat':
                            sdir=sdir[:-1]+'_filtered/'    
                    if not os.path.isdir(sdir):
                        print(f"make dir: {sdir}")
                        os.makedirs(f"{sdir}")
                    
                    print('calculating measurement covariance matrix')
                    measurement_covariance=[]
                    for i, dataset in enumerate(datasets):
                        if dataset.error_inflation:
                            measurement_covariance.extend(dataset.get_measurement_error_inflation_array(measurement_error[i],data[i]))
                        else:
                            measurement_covariance.extend([measurement_error[i]]*dataset.size)
                    measurement_covariance=np.array(measurement_covariance)
                    print('calculated measurement covariance matrix')
                    # run inverion, save dataset
                    spath=f"{sdir}/{start_date.strftime('%Y%m%d')}-{end_date.strftime('%Y%m%d')}_{bg_ds}_bg.nc"
                    if os.path.isfile(spath):
                        print('file already exists')
                        print(f'check {spath}')
                        break
                    if WITH_FLAT:
                        run_inv_TM5_prior_flat(flat_prior, TM5_prior, prior_cov, measurements,measurement_covariance, footprints,spath,SAVE_AK, footprint_col_name=f_col)
                    else:
                        run_inv(TM5_prior, prior_cov, measurements,measurement_covariance, footprints,spath,SAVE_AK, footprint_col_name=f_col)
                    del measurement_covariance
                    if FILTER_GOSAT_MEAS:
                        # save with pointspec_to_keep
                        merged.to_csv(f'{sdir}/gosat_meas_pointspec_to_keep.nc')                        
                        
                del footprints
            del prior_cov
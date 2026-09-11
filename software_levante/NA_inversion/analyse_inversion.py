import xarray as xr
import os
import pandas as pd
import geopandas
from shapely.geometry import Polygon,box
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from utils import get_start_date_of_week, get_unique_time
from itertools import product


class measurement_dataset():
    """Measurement dataset for which Flexpart was run

    Args:
        name: name of dataset
        version: version of dataset
        flexpart_path: Flexpart output directory
        footprints_path: path to footprint dataset
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
            footprints_path: str,
            num_parts: int, 
            measured: str,
            measurement_error,
            height_levels: str,
            pressure_weighted: bool = False, 
            averaging_kernel: bool = False
            ):
        
        self.name, self.version = name, version
        self.flexpart_path,self.footprints_path = flexpart_path,footprints_path
        self.num_parts=num_parts
        self.measurement_error = measurement_error
        self.measured, self.height_levels  = measured, height_levels
        self.pressure_weighted=pressure_weighted
        self.averaging_kernel=averaging_kernel



# postprocessing of inversion data
def get_co2_val_from_fluxes(start_date, end_date, dataset, inv_flux_path,BG, BG_molefrac_bg_str, res,ref_flux_dir, s_dir, f_str):
    ''' get co2 values from fluxes and footprints for combined GOSAT+IS inversion for fixed meas_err and land_err values, 
        using effektive footprints from f_str
        using unscaled footprints for reference fluxes
    Args:    
        start_date: start date of fluxes that are used
        end_date: end date of fluxes that are used
        dataset: measurement dataset
        inv_flux_path: path to inversion result fluxes
        BG_molefrac_bg_str: string defining the background, eg '{bg_str}', '395ppm' for fixed background value
        res: resolution, needed for reference fluxes
        ref_flux_dir: path to directory containing the reference fluxes
        s_dir: path to directory where output should be saved
        f_str: string defining which footprint to use
    '''
    # read footprint data
    print(dataset.footprints_path)
    footprints=xr.open_dataset(dataset.footprints_path)

    # reset pointspec dimension
    footprints['pointspec']=np.arange(0,footprints.pointspec.size)

    # read flux data
    flux_data=xr.open_dataset(inv_flux_path)
    ref_RemoteC_IS_flux=xr.open_dataset(f'{ref_flux_dir}/flux_{res}x{res}_RemoTeC_2.4.1+IS_cut_weekly.nc').sel(time=slice(start_date+relativedelta(months=1), end_date+relativedelta(months=-1)))
    
    
    # calculate co2/xco2 values, save dataframes
    # value from flux*footprint, sum    
    footprints['prior_val']=(flux_data.TM5_prior_flux*footprints[f_str]).sum(dim=['time', 'latitude', 'longitude'])
    footprints['posterior_val']=(flux_data.TM5_posterior_flux*footprints[f_str]).sum(dim=['time', 'latitude', 'longitude'])

    # from TM5-4DVar fluxes, those with unscaled fluxes
    footprints['ref_RemoteC_IS_val']=(ref_RemoteC_IS_flux.total_flux*footprints.spec001_mr).sum(dim=['time', 'latitude', 'longitude'])
    
    footprints['prior_diff']=footprints['prior_val']-footprints[dataset.measured]+footprints[f'{BG}_{BG_molefrac_bg_str}_background']
    footprints['posterior_diff']=footprints['posterior_val']-footprints[dataset.measured]+footprints[f'{BG}_{BG_molefrac_bg_str}_background']
    footprints['ref_RemoteC_IS_diff']=footprints['ref_RemoteC_IS_val']-footprints[dataset.measured]+footprints[f'{BG}_{BG_molefrac_bg_str}_background']
    footprints=footprints.drop_vars(f_str)
    # save dataframes
    footprints.to_netcdf(f'{s_dir}/{dataset.name}_data_{res}x{res}_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{BG_molefrac_bg_str}_bg.nc')
    print(f'saved files to {s_dir}')

def get_high_res_TM5_co2(month_date, prior_flux_path, f_dir, s_dir, cols, col_name='prior_val'):
    ''' calculate high resolution co2 values from TM5-4DVar 2x3° flux and high resolution footprints (1x1° and hourly)
    Args:
        month_date: date for the selected month
        prior_flux_path: path to high res prior
        f_dir: path to footprint directory
        s_dir: path to output directory
        cols: list of data columns from footprint to add to saved dataset
        col_name: data column name, defaults to 'prior_val'
    Returns:
        nothing, saves values as dataframe'''
    # TM5-4DVar prior
    prior_flux=xr.open_dataset(prior_flux_path)[['grid_cell_area', 'total_flux']]
    prior_flux=prior_flux.sel(time=slice(month_date-relativedelta(months=1), month_date))
    # high resolution footprints
    footprints=xr.open_dataset(f'{f_dir}/high_res_scaled_footprints_{month_date.strftime("%Y_%m")}.nc')
    # reindex prior to same time
    prior_flux=prior_flux.reindex(time=footprints.time, method='ffill')
    # get value
    val=(prior_flux['total_flux']*footprints['spec001_mr']).to_dataset(name=col_name)
    # sum over time and space
    val=val.sum(dim=['time', 'latitude', 'longitude'])
    # add meas data 
    for col in cols:
        val[col]=footprints[col]
    spath=f'{s_dir}/high_res_prior_{month_date.strftime("%Y_%m")}.nc'
    val.to_netcdf(spath)
    print(f'saved to: {spath}')
    return

def get_weekly_meas_range(data_dir, start_date, end_date, BG, bg_str):
    ''' Get range of measurements and their standard deviation during a week, for insitu, GOSAT measurements and modeled values
    Args:
        data_dir: path to inversion output directory
    Returns:
        nothing, saves weekly mean, std and range of measurement data for each grid cell
    '''
    # read meas data
    gosat_weekly = xr.open_dataset(f'{data_dir}/gosat_data_2x2_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc')
    is_weekly = xr.open_dataset(f'{data_dir}/insitu_data_2x2_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc')

    #  get enhancement + bg
    is_weekly['prior_val+bg']=is_weekly.prior_val+is_weekly[f'{BG}_{bg_str}_background']
    is_weekly['posterior_val+bg']=is_weekly.posterior_val+is_weekly[f'{BG}_{bg_str}_background']
    gosat_weekly['prior_val+bg']=gosat_weekly.prior_val+gosat_weekly[f'{BG}_{bg_str}_background']
    gosat_weekly['posterior_val+bg']=gosat_weekly.posterior_val+gosat_weekly[f'{BG}_{bg_str}_background']

    # get release week and release box
    lat_list=np.arange(13,56,step=2.)
    lon_list=np.arange(-133,-62,step=2.)
    
    gosat_weekly['release_week']=(('pointspec'), [get_start_date_of_week(pd.to_datetime(gosat_weekly['release_time'][i].item()).replace(hour=0, minute=0, second=0)) for i in range(gosat_weekly.pointspec.size)])
    gosat_weekly['release_box_lat']=(('pointspec', ), [lat_list[np.abs(gosat_weekly.release_lat.values[i]-lat_list).argmin()].item() for i in range(gosat_weekly.pointspec.size)])
    gosat_weekly['release_box_lon']=(('pointspec', ), [lon_list[np.abs(gosat_weekly.release_lon.values[i]-lon_list).argmin()].item() for i in range(gosat_weekly.pointspec.size)])
    is_weekly['release_week']=(('pointspec'), [get_start_date_of_week(pd.to_datetime(is_weekly['release_time'][i].item()).replace(hour=0, minute=0, second=0)) for i in range(is_weekly.pointspec.size)])
    is_weekly['release_box_lat']=(('pointspec', ), [lat_list[np.abs(is_weekly.release_lat.values[i]-lat_list).argmin()].item() for i in range(is_weekly.pointspec.size)])
    is_weekly['release_box_lon']=(('pointspec', ), [lon_list[np.abs(is_weekly.release_lon.values[i]-lon_list).argmin()].item() for i in range(is_weekly.pointspec.size)])
    # assign coordinates
    gosat_weekly=gosat_weekly.assign_coords(release_week=gosat_weekly.release_week, release_box_lat=gosat_weekly.release_box_lat, release_box_lon=gosat_weekly.release_box_lon)
    is_weekly=is_weekly.assign_coords(release_week=is_weekly.release_week, release_box_lat=is_weekly.release_box_lat, release_box_lon=is_weekly.release_box_lon)

    # groupby week, get mean and std
    # mean
    gosat_weekly_mean=gosat_weekly.groupby(['release_week', 'release_box_lat', 'release_box_lon']).mean()
    is_weekly_mean=is_weekly.groupby(['release_week', 'release_box_lat', 'release_box_lon']).mean()    
    
    # columns to get range and std for
    gosat_cols=[f'{BG}_{bg_str}_xco2_interpolated', 'xco2','prior_val','prior_diff', 'prior_val+bg','posterior_val', 'posterior_diff', 'posterior_val+bg', f'{BG}_{bg_str}_background', f'{BG}_molefrac_diff_interpolated']
    is_cols=[f'{BG}_{bg_str}_co2_interpolated', 'co2_val[ppm]','prior_val','prior_diff', 'prior_val+bg','posterior_val','posterior_diff',  'posterior_val+bg',f'{BG}_{bg_str}_background', f'{BG}_molefrac_diff_interpolated']
    for col in gosat_cols:
        gosat_weekly_mean[f'{col}_min']=gosat_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).min()
        gosat_weekly_mean[f'{col}_max']=gosat_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).max()
        gosat_weekly_mean[f'{col}_std']=gosat_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).std()
    for col in is_cols:
        is_weekly_mean[f'{col}_min']=is_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).min()
        is_weekly_mean[f'{col}_max']=is_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).max()
        is_weekly_mean[f'{col}_std']=is_weekly[col].groupby(['release_week', 'release_box_lat', 'release_box_lon']).std()
    # select measurements in 2010
    gosat_weekly_mean=gosat_weekly_mean.sel(release_week=slice(dt.date(2021,1,1), dt.date(2021,12,31)))
    is_weekly_mean=is_weekly_mean.sel(release_week=slice(dt.date(2021,1,1), dt.date(2021,12,31)))
    # save
    if not os.path.isdir(f'{data_dir}/weekly_data/'):
        os.makedirs(f'{data_dir}/weekly_data/')
    gosat_weekly_mean.to_netcdf(f'{data_dir}/weekly_data/gosat_data_weekly_mean_2010_RemoTeC_2.4.1+IS_bg.nc')
    is_weekly_mean.to_netcdf(f'{data_dir}/weekly_data/is_data_weekly_mean_2010_RemoTeC_2.4.1+IS_bg.nc')
    print(f'saved weekly meas data for directory: {data_dir}')
    return

# plots
def plot_co2_xco2_diff_hist(gosat_path, is_path, s_dir,bg_str,cols=['posterior'], title_str=''):
    ''' plot difference between modeled values and measurements for insitu and gosat measurements
    Args:
        gosat_path: path to gosat xco2 values
        is_path: path to insitu co2 values
        s_dir: path to directory where the plots should be saved
        bg_str: string indicating dataset used for the background
        cols: list of columns that should be plotted, from ['prior','posterior'], defaults to posterior
        title_str: title string for the plot, defaults to empty string
    Returns:
        nothing, saves plots to s_dir
    '''
    # read data
    gosat_data=xr.open_dataset(gosat_path)
    insitu_data=xr.open_dataset(is_path)
    # get max val from 'TM5_molefrac_diff', 'prior_diff', 'posterior_diff'
    gosat_max=np.ceil(np.max(np.abs(gosat_data[['ref_RemoteC_IS_diff', 'prior_diff', 'posterior_diff']]).to_array()).item())
    insitu_max=np.ceil(np.max(np.abs(insitu_data[['ref_RemoteC_IS_diff', 'prior_diff', 'posterior_diff']]).to_array()).item())
    # def bins
    bin_size=1
    gosat_bins=np.arange(-gosat_max, gosat_max, step=bin_size/2)
    insitu_bins=np.arange(-insitu_max, insitu_max, step=bin_size)
    # plot
    fig, ax=plt.subplots(1,2,figsize=(12,5))
    fig.suptitle(title_str)
    ax[0].hist(insitu_data.ref_RemoteC_IS_diff, insitu_bins,color='lightgreen',  label=fr'TM5-4DVar, $\overline{{d}}={np.round(np.mean(insitu_data.CAMS_molefrac_diff_interpolated).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(insitu_data.CAMS_molefrac_diff_interpolated)).compute().item(),2)}$')
    ax[0].set_xlabel('co2 difference to measurements [ppm]')
    ax[0].set_title(f'In-situ measurements')
    ax[1].hist(gosat_data.ref_RemoteC_IS_diff, gosat_bins,color='lightgreen',  label=fr'TM5-4DVar, $\overline{{d}}={np.round(np.mean(gosat_data.CAMS_molefrac_diff_interpolated).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(gosat_data.CAMS_molefrac_diff_interpolated)).compute().item(),2)}$')
    ax[1].set_xlabel('xco2 difference to measurements [ppm]')
    ax[1].set_title(f'Gosat RemoTeC v2.4.1 measurements')
    if 'prior' in cols:
        ax[0].hist(insitu_data.prior_diff, insitu_bins, color='orange', alpha=0.5,  label=fr'prior, $\overline{{d}}={np.round(np.mean(insitu_data.prior_diff).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(insitu_data.prior_diff)).compute().item(),2)}$')
        ax[1].hist(gosat_data.prior_diff, gosat_bins, color='orange', alpha=0.5,  label=fr'prior, $\overline{{d}}={np.round(np.mean(gosat_data.prior_diff).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(gosat_data.prior_diff)).compute().item(),2)}$')
    if 'posterior' in cols:
        ax[0].hist(insitu_data.posterior_diff, insitu_bins, color='red', alpha=0.5,  label=fr'posterior, $\overline{{d}}={np.round(np.mean(insitu_data.posterior_diff).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(insitu_data.posterior_diff)).compute().item(),2)}$')
        ax[1].hist(gosat_data.posterior_diff, gosat_bins, color='red', alpha=0.5,  label=fr'posterior, $\overline{{d}}={np.round(np.mean(gosat_data.posterior_diff).compute().item(),2)}, \overline{{|d|}}={np.round(np.mean(np.abs(gosat_data.posterior_diff)).compute().item(),2)}$')
    for i in range(0,len(ax)):
        ax[i].grid(alpha=0.5)
        ax[i].legend(fontsize='small')
    ax[0].set_xlim(-30,30)
    ax[1].set_xlim(-15,15)
    fig.tight_layout()
    sfile_name='hist_diff_comp'
    for c in cols:
        sfile_name+=f'_{c}'
    print(f'saving {s_dir}/{sfile_name}_{bg_str}.png')
    fig.savefig(f'{s_dir}/{sfile_name}_{bg_str}.png')
    plt.close()

def plot_4x4_flux_overwiev(data_dir,start_date, end_date,ref_flux_dir, bg_ds,title_str, plot_mean_2x2=False, plot_entrie_period=True):
    ''' for 4x4° inversion, plot flux timeseries
    Args:
        data_dir: path to inversion flux data with correlation
        start_date, end_date: dt.datetime defing time period
        ref_flux_dir: TM5-4DVar flux directory
        bg_ds: dataset used for the background calculation
        title_str: string used as suptitle for the plot
        plot_mean_2x2: set to True if should plot the mean of the respective 2x2 fluxes
        plot_entrie_period: set to False to only plot fluxes for 2010
    Returns:
        nothing, saves plot
    '''
    
    res= 4
    interp_region=[20,48,-126,-70]       # lat_min, lat_max, lon_min, lon_max
    # TM5-4DVar reference fluxes
    RemoTeC_IS_ref_data=xr.open_dataset(f'{ref_flux_dir}/flux_4x4_RemoTeC_2.4.1+IS_cut_weekly.nc').sel(time=slice(start_date+relativedelta(months=-1), end_date+relativedelta(months=1)))

    # path to 4x4 flux data
    inv_flux_path= f'{data_dir}/{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_ds}_bg.nc'
    data_4x4=xr.open_dataset(inv_flux_path)

    # 4x4 flux data without correlation
    data_4x4_without=xr.open_dataset(inv_flux_path.replace('with_correlation', 'no_correlation'))

    if plot_mean_2x2:
        # 2x2 inversions
        inv_flux_path_2x2=inv_flux_path.replace('4x4', '2x2')
        data_2x2=xr.open_dataset(inv_flux_path_2x2)

    lat_lon_list=[(lat, lon) for lat in data_4x4.latitude.values for lon in data_4x4.longitude.values]
    # sort by latitude descending, longitude ascending
    lat_lon_list=sorted(lat_lon_list, key=lambda x: (-x[0], x[1]))
    # plot overview
    fig, ax=plt.subplots(data_4x4.latitude.size, data_4x4.longitude.size, figsize=(data_4x4.longitude.size*6, data_4x4.latitude.size*6), sharey=True)
    ax=ax.flatten()
    # marker size
    ms=2
    for i in range(0, len(lat_lon_list)):
        lat_sel, lon_sel =lat_lon_list[i]
        # select region
        RemoTeC_IS_ref_data_sel=RemoTeC_IS_ref_data.sel(latitude=lat_sel, longitude=lon_sel)
        data_4x4_sel=data_4x4.sel(latitude=lat_sel, longitude=lon_sel)
        data_4x4_without_sel=data_4x4_without.sel(latitude=lat_sel, longitude=lon_sel)
        if plot_mean_2x2:
            data_2x2_sel=data_2x2.sel(latitude=slice(lat_sel-2, lat_sel+2), longitude=slice(lon_sel-2,lon_sel+2)).mean(dim=['latitude','longitude'])
        
        # plot prior and posterior
        ax[i].plot(data_4x4_sel.time, data_4x4_sel.TM5_prior_flux, c='orange',marker='o',markersize=ms, label='TM5-4DVar prior')
        ax[i].fill_between(data_4x4_sel.time, data_4x4_sel.TM5_prior_flux-data_4x4_sel.prior_uncertainty, data_4x4_sel.TM5_prior_flux+data_4x4_sel.prior_uncertainty, color='orange', alpha=0.2)
        ax[i].plot(data_4x4_sel.time, data_4x4_sel.TM5_posterior_flux, c='r', marker='o', markersize=ms, label='4x4 posterior')
        ax[i].fill_between(data_4x4_sel.time, data_4x4_sel.TM5_posterior_flux-data_4x4_sel.TM5_posterior_std, data_4x4_sel.TM5_posterior_flux+data_4x4_sel.TM5_posterior_std, color='r', alpha=0.2)
        if not plot_mean_2x2:
            ax[i].plot(data_4x4_without_sel.time, data_4x4_without_sel.TM5_posterior_flux, c='brown', marker='o', markersize=ms, label='4x4 posterior, no correlation')
            ax[i].fill_between(data_4x4_without_sel.time, data_4x4_without_sel.TM5_posterior_flux-data_4x4_without_sel.TM5_posterior_std, data_4x4_without_sel.TM5_posterior_flux+data_4x4_without_sel.TM5_posterior_std, color='brown', alpha=0.2)
        if plot_mean_2x2:
            ax[i].plot(data_2x2_sel.time, data_2x2_sel.TM5_posterior_flux, c='green', marker='o', markersize=ms, label='mean 2x2 posterior')
        # plot reference
        ax[i].plot(RemoTeC_IS_ref_data_sel.time, RemoTeC_IS_ref_data_sel.total_flux,c='b', marker='x', markersize=ms, label='reference: RemoTeC_2.4.1+IS')
        
        ax[i].legend()
        ax[i].set_title(f'lat/lon: {lat_sel}/{lon_sel}')
        if plot_entrie_period:
            ax[i].set_xlim(start_date, end_date)
        else:
            ax[i].set_xlim(dt.date(2010,1,1), dt.date(2010,12,31))
        ax[i].tick_params(labelleft=True)
        ax[i].set_ylim(-2.5e-7, 2e-7)
        
        # You can change the step of range() as you prefer (now, it selects each third month) 
        ax[i].xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,3)))
        ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=range(1,13,1)))
        ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax[i].grid(which='both', c='grey', alpha=0.5)
    # box around interpretation region
    rect = plt.Rectangle(
        # (lower-left corner), width, height
        (1/9, 2/11*0.98), 7/9, 7/11*0.98, fill=False, color="lightgreen", lw=6,# alpha=0.5, 
        zorder=1000, transform=fig.transFigure, figure=fig
    )
    fig.patches.extend([rect])

    fig.suptitle(title_str)
    fig.tight_layout(rect=[0, 0., 1, 0.98])
    if plot_mean_2x2:
        spath_fig=f'{data_dir}/flux_comp_4x4_mean2x2_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_ds}.png'
    else:
        spath_fig=f'{data_dir}/flux_comp_4x4_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_ds}.png'
    if not plot_entrie_period:
        spath_fig=spath_fig.replace('.png', '_only2010.png')
    fig.savefig(spath_fig)
    print(f'saved to: {spath_fig}')
    plt.close()

def plot_inner_2x2_region_fluxes(data_dir,ref_flux_dir, gosat_err_val, meas_err_val,corr_str, color_sel,start_date, end_date, bg_str = 'RemoTeC_2.4.0+IS'):
    ''' plot time series of prior, posterior and TM5-4DVar reference fluxes for all grid cells in inner domain
    Args:
        data_dir:   path to inversion directory
        ref_flux_dir: TM5-4DVar flux directory
        gosat_err_val: gosat measurement error value
        meas_err_val: insitu measurement error value
        corr_str: sting indicating if with or withour correlation
        color_sel: color for posterior flux
        bg_str: TM5-4DVar dataset used for the background, defaults to RemoTeC_2.4.0+IS
    Returns: nothing, saves plot'''
    
    # read data
    data=xr.open_dataset(f'{data_dir}/{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc').sel(latitude=slice(20,48), longitude=slice(-123,-73))
    # read reference
    bg_str2 = 'RemoTeC_2.4.1+IS'
    ref_flux_2x2=xr.open_dataset(f'{ref_flux_dir}/flux_2x2_{bg_str2}_cut.nc').sel(latitude=slice(20,48), longitude=slice(-123,-73))

    title_str=f'Total 2x2 flux, {corr_str} correlation, with total scaling,{gosat_err_val}pp gosat meas err, {meas_err_val}ppm measurement error'
    lat_lon_list=[(lat, lon) for lat in data.latitude.values for lon in data.longitude.values]
    # sort by latitude descending, longitude ascending
    lat_lon_list=sorted(lat_lon_list, key=lambda x: (-x[0], x[1]))
    # plot overview
    fig, ax=plt.subplots(data.latitude.size, data.longitude.size, figsize=(data.longitude.size*6, data.latitude.size*6), sharey=True)
    ax=ax.flatten()
    # marker size
    ms=2
    for i in range(0, len(lat_lon_list)):
        lat_sel, lon_sel =lat_lon_list[i]
        # select region
        ref_flux_2x2_sel=ref_flux_2x2.sel(latitude=lat_sel, longitude=lon_sel)
        data_sel=data.sel(latitude=lat_sel, longitude=lon_sel)
        # plot prior
        ax[i].plot(data_sel.time, data_sel.TM5_prior_flux, c='orange', marker='o', markersize=ms, label='prior')
        ax[i].fill_between(data_sel.time, data_sel.TM5_prior_flux-data_sel.prior_uncertainty, data_sel.TM5_prior_flux+data_sel.prior_uncertainty, color='orange', alpha=0.2)        
        # plot posterior
        ax[i].plot(data_sel.time, data_sel.TM5_posterior_flux, c=color_sel, marker='o', markersize=ms, label='posterior')
        ax[i].fill_between(data_sel.time, data_sel.TM5_posterior_flux-data_sel.TM5_posterior_std, data_sel.TM5_posterior_flux+data_sel.TM5_posterior_std, color=color_sel, alpha=0.2)        
        # plot reference flux
        ax[i].plot(ref_flux_2x2_sel.time, ref_flux_2x2_sel.total_flux,c='b', marker='x', markersize=ms, label='reference: RemoTeC_2.4.1+IS')
        
        ax[i].legend()
        ax[i].set_title(f'lat/lon: {lat_sel}/{lon_sel}')
        ax[i].set_xlim(start_date,end_date)
        #ax[i].set_xlim(dt.date(2010,1,1), dt.date(2010,12,31))
        ax[i].tick_params(labelleft=True)
        ax[i].set_ylim(-4e-7, 4e-7)
        
        # You can change the step of range() as you prefer (now, it selects each third month) 
        ax[i].xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,3)))
        ax[i].xaxis.set_minor_locator(mdates.MonthLocator(bymonth=range(1,13,1)))
        ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax[i].grid(which='both', c='grey', alpha=0.5)
    fig.suptitle(title_str)
    fig.tight_layout(rect=[0, 0., 1, 0.98])
    spath_fig=f'{data_dir}/flux_comp_2x2_inner_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}.png'
    fig.savefig(spath_fig)
    print(f'saved to: {spath_fig}')
    plt.close()
    return


# main
if __name__ == "__main__":
    # TODO adapt paths and variables before running!
    GET_EAST_WEST_MASKS=False
    GET_CO2_VALS=True              # calculate modeled measurement values from our posterior
    GET_HIGH_RES_PRIOR_CO2=False            # calculate modeled measurement values from TM5-4DVar prior and high resolution flexpart footprints (1x hourly)
    GET_HIGH_RES_TM5_POSTERIOR_CO2=False    # same for TM5-4DVar posterior flux
    GET_WEEKLY_MEAS=False           # get weekly mean and range of measurements/ modeled values
    
    PLOT_DIFF_HIST=False
    PLOT_INNER_2X2_FLUXES=False     # plot fluxes for 2x2 inversion, inner domain
    PLOT_4x4_FLUX_MAP=False     # plot fluxes for 4x4 inversion
    
    start_date = dt.date(2020,10,1)
    end_date = dt.date(2022,3,31)
    
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

    # path to Flexpart directory with subdirectories for datasets
    flex_dir=directory+'/Flexpart/' 
    # footprint variable name
    f_str='spec001_mr_scaled'    # 'spec001_mr', 'spec001_mr_scaled'
    footprint_subdirectory= "prep_footprints_TM5_3_years/scaled_weekly/" 
    inversion_subdirectory="inversions_new/TM5_diurnal/scaled/"
    # with/without correlation
    corr_str_list=['with'] #['no', 'temporal', 'spatial', 'temporal_1M', 'with_1M', 'with_e04']          # 'no', 'with', 'temporal', 'spatial', 'temporal_1M', 'with_1M', 'with_e04'
    # spatial resolution
    res=2
    
    # list of errors
    gosat_err_list=[2]          # [0.25,0.5, 0.8,1, 1.2, 1.5,2,2.5,3]
    is_err_list=[3,5]             # [0.5,1,2,4,6,8,10]
    # selected error values
    gosat_err_val_sel=2
    is_err_val_sel=3

    
    #path to reference fluxes (TM5-4DVar flux directory)
    ref_flux_dir = directory+'/TM54DVar/TM54DVar_fluxes_3years/'   

    # measurement datasets:
    gosat=measurement_dataset(
        name = 'gosat', 
        version = '2.4.1', 
        flexpart_path = flex_dir+'/'+'RemoTeCv240',
        footprints_path = flex_dir+'/'+'RemoTeCv240'+'/'+footprint_subdirectory+f'/high_res_scaled_footprints_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_2x2_weekly.nc',
        num_parts = 40000,
        measured = 'xco2',
        measurement_error=gosat_err_list,
        height_levels = 'continuous profile',
        pressure_weighted = True,
        averaging_kernel = True)

    insitu=measurement_dataset(
        name='insitu',
        version='GLOBALVIEWplus_v10.1',
        flexpart_path=flex_dir+'insitu',
        footprints_path = flex_dir+'/'+'insitu'+'/'+footprint_subdirectory+f'/high_res_scaled_footprints_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_2x2_weekly.nc',
        num_parts=40000,
        measured = 'co2', # old: 'co2_val[ppm]',
        measurement_error=is_err_list,
        height_levels = 'single point')

    tccon=measurement_dataset(
        name='TCCON',
        version='',
        flexpart_path=flex_dir+'TCCON',
        footprints_path = flex_dir+'/'+'TCCON'+'/'+footprint_subdirectory+f'/high_res_scaled_footprints_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_2x2_weekly.nc',
        num_parts=40000,
        measured = 'xco2',
        height_levels = 'multiple points',
        measurement_error='',
        pressure_weighted=True,
        averaging_kernel=True)
    
    datasets_in_inversion=[gosat,insitu]
    
    measurement_errors=[]
    for i, dataset in enumerate(datasets_in_inversion):
        measurement_errors.append(dataset.measurement_error)
    measurement_errors=list(product(*measurement_errors))

    if GET_CO2_VALS:
        for dataset in [tccon,insitu, gosat]:
            print(dataset.name)
            for corr_str in corr_str_list:
                print(corr_str)
                for measurement_error in measurement_errors:
                    # create directory in which footprints will be saved
                    flux_dir=f'{flex_dir+inversion_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_str}/'
                    for i, inv_dataset in enumerate(datasets_in_inversion):
                        flux_dir+=f'{measurement_error[i]}ppm_{inv_dataset.name}_meas_err/'
                    flux_path=f'{flux_dir}{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{BG_molefrac_bg_str}_bg.nc'
                    get_co2_val_from_fluxes(start_date, end_date, dataset, flux_path,BG, BG_molefrac_bg_str, res,ref_flux_dir, flux_dir, f_str)

    if GET_HIGH_RES_PRIOR_CO2:
        # get high resolution prior
        prior_flux_path=f'{ref_flux_dir}/flux_1x1_prior_cut.nc'
        # insitu
        f_dir=f'{parent_data_dir}/insitu/prep_footprints/scaled_hourly'
        s_dir= f'{parent_data_dir}/insitu/prep_footprints/prior_co2_val'
        cols=['release_time', 'release_lat', 'release_lon','co2_val[ppm]',f'{BG}_{bg_str}_background', f'{BG}_{bg_str}_co2']
        if not os.path.isdir(s_dir):
            os.makedirs(s_dir)
        for date in pd.date_range(start_date, end_date, freq='MS'):
            get_high_res_TM5_co2(date, prior_flux_path, f_dir, s_dir, cols)
        # for gosat
        f_dir=f'{parent_data_dir}/RemoTeCv240/prep_footprints/scaled_hourly'
        s_dir= f'{parent_data_dir}/RemoTeCv240/prep_footprints/prior_co2_val'
        if not os.path.isdir(s_dir):
            os.makedirs(s_dir)
        cols=['release_time', 'release_lat', 'release_lon','xco2', 'xco2_err',f'{BG}_{bg_str}_background', f'{BG}_{bg_str}_xco2']
        for date in pd.date_range(start_date, end_date, freq='MS'):
            get_high_res_TM5_co2(date, prior_flux_path, f_dir, s_dir, cols)
    if GET_HIGH_RES_TM5_POSTERIOR_CO2:
        # same for posterior
        posterior_flux_path=f'{ref_flux_dir}/flux_1x1_RemoTeC_2.4.1+IS_cut.nc'
        # insitu
        f_dir=f'{parent_data_dir}/insitu/prep_footprints/scaled_hourly'
        s_dir= f'{parent_data_dir}/insitu/prep_footprints/posterior_RemoTeC_2.4.0+IS_co2_val'
        cols=['release_time', 'release_lat', 'release_lon','co2_val[ppm]',f'{BG}_{bg_str}_background', f'{BG}_{bg_str}_co2']
        if not os.path.isdir(s_dir):
            os.makedirs(s_dir)
        for date in pd.date_range(start_date, end_date, freq='MS'):
            get_high_res_TM5_co2(date, posterior_flux_path, f_dir, s_dir, cols, col_name='RemoTeC_2.4.0+IS_val')
        # for gosat
        f_dir=f'{parent_data_dir}/RemoTeCv240/prep_footprints/scaled_hourly'
        s_dir= f'{parent_data_dir}/RemoTeCv240/prep_footprints/posterior_RemoTeC_2.4.0+IS_co2_val'
        if not os.path.isdir(s_dir):
            os.makedirs(s_dir)
        cols=['release_time', 'release_lat', 'release_lon','xco2', 'xco2_err',f'{BG}_{bg_str}_background', f'{BG}_{bg_str}_xco2']
        for date in pd.date_range(start_date, end_date, freq='MS'):
            get_high_res_TM5_co2(date, posterior_flux_path, f_dir, s_dir, cols, col_name='RemoTeC_2.4.0+IS_val')
    
    if GET_WEEKLY_MEAS:
        for corr_str in corr_str_list:
            # fixed gosat_err_val, variable meas_err_val
            gosat_err_val=gosat_err_val_sel
            for meas_err_val in is_err_list:      #
                inv_dir=f'{parent_data_dir}/{inversions_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
                get_weekly_meas_range(inv_dir, start_date, end_date,BG,bg_str)
            meas_err_val=is_err_val_sel
            for gosat_err_val in gosat_err_list:        #
                inv_dir=f'{parent_data_dir}/{inversions_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
                get_weekly_meas_range(inv_dir, start_date, end_date,BG,bg_str)

    if PLOT_INNER_2X2_FLUXES:
        gosat_err_val=gosat_err_val_sel
        cmap = mpl.colormaps['Reds']
        # Take colors at regular intervals spanning the colormap.
        colors = cmap(np.linspace(0.5, 0.9, len(is_err_list)))
        for i in range(len(is_err_list)):
            meas_err_val=is_err_list[i]
            inv_dir=f'{parent_data_dir}/{inversions_subdirectory}/2x2/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
            plot_inner_2x2_region_fluxes(inv_dir,ref_flux_dir, gosat_err_val, meas_err_val,corr_str,colors[i],start_date, end_date, bg_str = bg_str)#'RemoTeC_2.4.0+IS')
        
        meas_err_val=is_err_val_sel
        cmap = mpl.colormaps['Reds']
        # Take colors at regular intervals spanning the colormap.
        colors = cmap(np.linspace(0.5, 0.9, len(gosat_err_list)))
        for i in range(len(gosat_err_list)):
            gosat_err_val=gosat_err_list[i]
            inv_dir=f'{parent_data_dir}/{inversions_subdirectory}/2x2/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
            plot_inner_2x2_region_fluxes(inv_dir,ref_flux_dir, gosat_err_val, meas_err_val,corr_str,colors[i], start_date, end_date, bg_str = bg_str) #'RemoTeC_2.4.0+IS')
    if PLOT_4x4_FLUX_MAP:        
        for meas_err_val in is_err_list:           
            data_dir=f'{parent_data_dir}/{inversions_subdirectory}/4x4/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
            title_str=f'Total CO2 flux, {meas_err_val}ppm insitu meas_err,{gosat_err_val}ppm gosat meas_err, {bg_str} background, {f_str}'
            # plot_4x4_flux_overwiev(data_dir,start_date, end_date,ref_flux_dir, bg_str,title_str)
            plot_4x4_flux_overwiev(data_dir,start_date, end_date,ref_flux_dir, bg_str,title_str, plot_entrie_period=False, plot_mean_2x2=True)
    if PLOT_DIFF_HIST:
        for corr_str in corr_str_list:  
            # fixed gosat_err_val, variable meas_err_val
            gosat_err_val=gosat_err_val_sel
            for meas_err_val in is_err_list:
                flux_dir=f'{parent_data_dir}/{inversions_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
                gosat_path=f'{flux_dir}/gosat_data_{res}x{res}_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc'
                insitu_path=f'{flux_dir}/insitu_data_{res}x{res}_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc'
                # plot difference histogramms
                plot_co2_xco2_diff_hist(gosat_path, insitu_path, flux_dir,bg_str,cols=['prior','posterior'], title_str=f'{res}x{res}, {corr_str} correlation, {gosat_err_val}ppm gosat meas_err, {meas_err_val}ppm insitu meas_err')
                # plot_co2_xco2_diff_hist(gosat_path, insitu_path, flux_dir,bg_str,cols=['prior'], title_str=f'{res}x{res}, {corr_str} correlation, {gosat_err_val}ppm gosat meas_err, {meas_err_val}ppm insitu meas_err')
            meas_err_val=is_err_val_sel 
            for gosat_err_val in gosat_err_list:
                flux_dir=f'{parent_data_dir}/{inversions_subdirectory}/{res}x{res}/{corr_str}_correlation/footprint_{f_str}/{gosat_err_val}ppm_gosat_meas_err/{meas_err_val}ppm_insitu_meas_err/'
                gosat_path=f'{flux_dir}/gosat_data_{res}x{res}_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc'
                insitu_path=f'{flux_dir}/insitu_data_{res}x{res}_{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}_{bg_str}_bg.nc'
                # plot difference histogramms
                plot_co2_xco2_diff_hist(gosat_path, insitu_path, flux_dir,bg_str,cols=['prior','posterior'], title_str=f'{res}x{res}, {corr_str} correlation, {gosat_err_val}ppm gosat meas_err, {meas_err_val}ppm insitu meas_err')




import xarray as xr
import os
import pandas as pd
import geopandas
import datetime as dt
import numpy as np
from shapely.geometry import Polygon


def read_grid_time_file(path):
    ''' Open FLEXPART grid_time file
    Args: 
        path: path to FLEXPART output directory
    Returns:
        dataset '''
    for file in os.listdir(path):
        if file.startswith('grid_time_'):
            filepath=path+file
            print('filepath: {}'.format(filepath))
    return xr.open_dataset(filepath, chunks="auto") 

def check_succesfull_footprints(start_date, end_date, data_dir):
    ''' Check if all footprints were calculated succesfully
    Args:
        start_date, end_date (datetime.date): start/end of releases that should be checked
        data_dir: path to directory containing the footprints with the subdirectories YYYY_mm/Release_YYYYmmdd/log.txt
    Returns:
        List of  dates for which the footprints were not calculated successfully
    '''
    date_list=[]
    for date in pd.date_range(start_date, end_date):
        dir_path=f'{data_dir}/{date.strftime("%Y_%m")}/Release_{date.strftime("%Y%m%d")}/'
        if os.path.isdir(dir_path):
            # dir_path
            with open(os.path.join(dir_path, "log.txt"), "r") as f:
                last_line = f.readlines()[-1]
                if not last_line==' CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLEXPART MODEL RUN!\n':
                    date_list.append(date)
                    print(f'Problem for Release_{date.strftime("%Y%m%d")}')
        else:
            print(f'no directory for {date.strftime("%Y%m%d")}')
    return date_list

def get_start_date_of_week(date):
    ''' Get monday date of the of the week containing the given date
    Args: date
    Returns: date of monday in corresponding week'''
    start_date=date-dt.timedelta(days=date.weekday())   # dt.date().weekday() goes from 0 to 6
    return start_date

def get_sim_start(release_dir):
    ''' Get start of simulation from COMMAND file in Release directory
    Args:
        release_dir: path to release directory
    Returns:
        sim_start as datetime object
    '''
    if not release_dir[-1]=='/':
        release_dir=release_dir+'/'
    
    f = open(f'{release_dir}COMMAND.namelist', "r")
    for line in f:
        if line.startswith(' IBDATE'):
            date_str=line[8:16]
        if line.startswith(' IBTIME'):
            time_str=line[8:14]    
    sim_start=dt.datetime.strptime(f'{date_str}_{time_str}', '%Y%m%d_%H%M%S')
    return sim_start

def getAreaOfGrid(lat_min,lat_max):
    """Get Dataframe with Area dependend on Latitude for a 1°x1° grid, between lat_min and lat_max"""
    # von Eva kopiert, aus CreateAndSaveGeoDataframeGFAS_Pernila.py
    # etwas angepasst, lat_min und max eingeführt
    AreaLat = []
    # Lat = range(895,-905,-10)
    Lat = range(int(lat_max*10),int(lat_min*10-10),-10)
    for i in Lat:
        geom = [Polygon(zip([100,100,101,101],[i/10-0.5,i/10+0.5,i/10+0.5,i/10-0.5]))]
        GData = geopandas.GeoDataFrame({'data':[1]}, geometry=geom)
        GData.crs = 'epsg:4326'
        # from lat/long grid to meter based coordinate system to get area in m^2
        GData = GData.to_crs({'proj':'cea'})
        AreaLat.append(GData.geometry.area[0])
    dfArea = pd.DataFrame({'Area':AreaLat,'Lat':np.array(Lat)/10})
    return(dfArea)

def add_to_gosat_sounding_data(data, columns, path, spath=''):
    """Add data to gosat sounding positions file, eg background, enhancements
    Args:
        data: list of data arrays that should be added to file
        columns: column names for dataarrays
        gosat_path: path to gosat sounding positions file that data should be added to
        gosat_spath: path to where gosat file, now including data, should be saved, defaults to gosat_path
    Returns:
        nothing, saves data into file
    """
    print(path)
    gosat_data=pd.read_csv(path)
    # add data
    for i in range(0,len(columns)):
        gosat_data[columns[i]]=data[i]
    # save to spath
    if spath=='':
        spath=path
    print(f'saving file: {spath}')
    gosat_data.to_csv(spath, index=False)
    return

def haversine(lat1, lon1, lat2, lon2):
    ''' Haversine function for geodesic distance, determines great circle distance between two points on a sphere
        used for spatial correlation in prior covariance matrix
    Args:
        lat1, lon1: latitude & longitude of first point
        lat2, lon2: latitude & longitude of second point
    Returns: geodesic distance (float)
    '''
    R = 6371  # Radius of Earth in km
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

def h(
    p,
    p0 = 1013.25,
    h0 = 0,
    Tb = 288.15,
    R = 8.3144598,
    g = 9.80665,
    M = 0.0289644,
    ):
    """Calculate height from pressure based on barrometric height formula as described here: https://en.wikipedia.org/wiki/Barometric_formula
    Args:
        p (float): pressure for height calculation [hPa]
        p0 (float, optional): reference pressure [hPa]. Defaults to 1013.25
        h0 (float, optional): height of reference [m]. Defaults to 0.
        Tb (float, optional): reference temperature [K]. Defaults to 288.15.
        R (float, optional): universal gas constant [J/(mol*K)]. Defaults to 8.3144598.
        g (float, optional): gravitational acceleration [m/s^2]. Defaults to 9.80665.
        M (float, optional): molar mass of Earth's air [kg/mol]. Defaults to 0.0289644.

    Returns:
        float: height of pressure compared to pressure p0 at height h0
    """ 
    height =  -R*Tb/(g*M) * np.log(p/p0) + h0 # [m]
    return height

def p(
    h,
    Tb = 288.15,
    h0 = 0,
    p0 = 1013.25,
    R = 8.3144598,
    g = 9.80665,
    M = 0.0289644,
    ):
    """Calculate factor of barrometric height formula as described here: https://en.wikipedia.org/wiki/Barometric_formula

    Args:
        h (fleat): height for factor calculation [m]
        Tb (float, optional): reference temperature [K]. Defaults to 288.15.
        h0 (float, optional): height of reference [m]. Defaults to 0.
        p0 (float, optional): reference pressure [hPa]. Defaults to 1013.25
        R (float, optional): universal gas constant [J/(mol*K)]. Defaults to 8.3144598.
        g (float, optional): gravitational acceleration [m/s^2]. Defaults to 9.80665.
        M (float, optional): molar mass of Earth's air [kg/mol]. Defaults to 0.0289644.

    Returns:
        float: float: pressure at height h in hPa
    """    
    pressure = p0 * np.exp(-g * M * (h - h0)/(R * Tb))
    return pressure

def mean_error(dx, axis=None):
    n = dx.shape[axis] if axis is not None else dx.size
    return np.sqrt(np.sum(dx**2, axis=axis) / n**2)


# need to drop time dependence of release_time variable for gosat data
def get_unique_time(values):
    valid = values[~pd.isna(values)]
    unique = pd.unique(valid)
    if len(unique) > 1:
        raise ValueError(f"Multiple unique non-NaT values found: {unique}")
    # Ensure return type is datetime64[ns], even for NaT
    return np.datetime64(unique[0], 'ns') # if len(unique) == 1 else np.datetime64('NaT', 'ns')
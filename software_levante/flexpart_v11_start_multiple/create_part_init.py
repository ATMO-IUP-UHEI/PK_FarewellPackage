import numpy as np
import xarray as xr
import yaml
import click

###################################################
# run this script e.g. with 
# python create_part_init.py --config_path configs/part_init_config.yaml
###################################################

@click.command()
@click.option("--config_path", required=True, help="path to config.yaml file")
def main(config_path):
    config = read_config(config_path)
    outpath,num_part,num_layers,num_releases, lat,lon,release_time,pmin,pmax,nspecies,species_id,species_mass,kindz =(config['outpath'], config['num_part'],config['num_layers'],config['num_releases'], config['lat'], config['lon'], config['release_time'], config['pmin'], config['pmax'], config['nspecies'], config['species_id'], config['species_mass'], config['kindz'])
    # create and save dataset
    create_part_init(outpath,num_part,num_layers,num_releases, lat,lon,release_time,pmin,pmax,nspecies,species_id,species_mass,kindz)

    
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

# TODO add other options for part_init than total column at one location
def create_part_init(outpath,num_part,num_layers,num_releases, lat,lon,release_time,pmin,pmax,nspecies,species_id,species_mass,kindz):
    """Creates part_init.nc file that can be used as initial condition for FLEXPART_v11 run
    Args:
        outpath: path to directory where file will be saved
        num_part: total number of particles released for one measurement
        num_layers: number of height layers per measurement
        num_releases: number of measurements to be written in part_ini.nc file
        lat,lon:    position of release
        release_time: Release time of each particle in seconds after simulation start (IBDATE/IBTIME for forward runs, IEDATE/IETIME for backward runs, set in COMMAND file)
        pmin,pmax:  minimum and maximum pressure, particles equally distributed in pressure
        nspecies: number of species that should be released
        species_id: id of species in SPECIES directory in FLEXPART options directory, CO2=41
        species_mass:   molecular mass of species in g/mol
        kindz:  integer, Reference level: 1=above ground, 2=above sea level, 3 for pressure in hPa

    Returns:
        nothing, creates part_init.nc file that can be used as initial condition for FLEXPART_v11 run
    """ 
    if not (num_part % num_layers)==0:
        print('ERROR: number of particles can not be evenly distributed on number of layers')
        return
    # adjust numpart so that number of particles per layer is equal
    part_per_layer=int(num_part/num_layers)
    # get particle ids, ranges from 0 to number of particles per measurement release
    particle_id=np.arange(0,num_part*num_releases,step=1)
    
    # particles equally dirtributed in pressure
    particles=np.linspace(pmax,pmin,num_part)   
    # get height of particles, same for each release
    particles_h=h(particles)
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
        longitude[i*num_part:(i+1)*num_part]=lon[i]
        latitude[i*num_part:(i+1)*num_part]=lat[i]
        time[i*num_part:(i+1)*num_part]=release_time[i]
        height[i*num_part:(i+1)*num_part]=particles_h
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
                    num_releases=(num_releases,'number of releases'))
    )

    # save dataset
    ds.to_netcdf(outpath+'part_ic.nc')
    print('Written file {}part_ic.nc'.format(outpath))
    return

def read_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
        return config

#########################################################################################

if __name__ == '__main__':
    main()
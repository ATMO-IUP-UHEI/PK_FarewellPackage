# flexpart_v11_start_multiple

## Necessary preparations for Flexpart runs: 
### Select measurements inside selected region and time period from measurement datasets
- create_GOSATpositions.py
  - run with: \
  *python create_GOSATpositions.py --gosat_config_path gosat_config_RemoTeC241.yaml*  \
  for RemoTeC v.2.4.0 use configs/gosat_config_RemoTeC240.yaml \
  for RemoTeC v.2.4.1 use configs/gosat_config_RemoTeC241.yaml 
  - **adapt before running**:
    - adapt config file, e.g. configs/gosat_config_RemoTeC240.yaml
        - GOSAT_dir
        - outdir
        - Lat_min, Lat_max, Long_min, Long_max:
        - startdate, enddate
  - automatically creates subdirectories: outdir/YYYY_MM/outfile_YYYYMMDD.csv,
  - creates a .nc file with all GOSAT measurements and their avearging kernels inside the selected region for each day
    - has option to Average measurements if close together 
    - files contain variables: 
      ```
      -time                   (sounding_dim)
      -latitude               (sounding_dim)
      -longitude              (sounding_dim)
      -xco2                   (sounding_dim)
      -xco2_err               (sounding_dim)
      -xco2_averaging_kernel  (sounding_dim, layer_dim)   normalized column averaging kernel XCO2
      -pressure_levels        (sounding_dim, level_dim)
      -co2_column_apriori     (sounding_dim)
      -co2_profile_apriori    (sounding_dim, layer_dim)
      ```


- create_ISpositions.py
  - run with python create_ISpositions.py --config_path configs/IS_config.yaml
  - **adapt before running**:
      - adapt config file, e.g. configs/IS_config.yaml
          - path to combined ObsPack dataset
          - outdir
          - startdate, enddate
  - automatically creates subdirectories: outdir/YYYY_MM/outfile_YYYYMMDD.csv,
  - creates a .csv file for all in-situ measurements in selected region for each day,  
    - insitu-measurements taken from on {YYYY}_lat{min}_{max}_lon{min}_{max}sel_mean_combined.nc file (created with selectObsPackData.py skript)
  - files contain columns: 
    ```
    -time     (at start of 4h mean time period) 
    -file     (ObsPack file name)
    -latitude
    -longitude
    -co2,
    -elevation[masl],
    -intake_height[magl]
    ```

- create_TCCONpositions.py
  - run with: *python create_TCCONpositions.py --TCCON_config_path configs/TCCON_config.yaml*
  - **adapt before running**:
      - adapt config file, e.g. configs/TCCON_config.yaml
          - data_dir      path to folder with individual measurement files for different TCCON stations
          - sites_path    path to sites.json file, (should be in download of TCCON data)
          - Lat_min, Lat_max, Long_min, Long_max
          - outdir
          - startdate, enddate
  - automatically creates subdirectories: outdir/YYYY_MM/outfile_YYYYMMDD.csv,
  - creates a .nc file with all TCCON measurements and their avearging kernels inside the selected region for each day
  - files contain variables: 
      ```
      - latitude               (sounding_dim)
      - longitude              (sounding_dim) 
      - xco2                   (sounding_dim) 
      - xco2_error             (sounding_dim) 
      - xco2_averaging_kernel  (sounding_dim, layer_dim) 
      - ak_pressure            (sounding_dim, layer_dim)
      - co2_profile_apriori    (sounding_dim, layer_dim)
      - co2_column_apriori     (sounding_dim) 
      - integration_operator   (sounding_dim, layer_dim) 
      - ak_altitude            (sounding_dim, layer_dim) 
      - prior_altitude         (sounding_dim, layer_dim) 
      - local_time             (sounding_dim) 
      - site_id                (sounding_dim) 
      - name                   (sounding_dim) 
      - location               (sounding_dim) 
      - time                   (sounding_dim) 
      ```


### Prepare FLEXPART runs: prepare necessary files and folders
**adapt before running**:
  - adapt configs/options_dummy/OUTGRID, output grid for FLEXPART footprints
  - optionally: Other options_dummy/COMMAND file options, adapt particle output variables options_dummy/PARTOPTIONS, ...

- prepare_GOSATruns.py
  - run with: *python prepare_GOSATruns.py --config_path configs/options_config_RemoTeC240.yaml*
  - **adapt before running**:
      - adapt 
        - config file: configs/options_config_RemoTeC240.yaml
          - options_dummy_path   path to '.../PK_FarewellPackage/software_levante/flexpart_v11_start_multiple/configs/options_dummy'
          - output_path 
          - input_paths         path to ERA5 data
          - available_paths     path to AVALILABLE file for ERA5 data
          - part_init_config_path
          - sim_length
          - RemoTeC_version
          - satellites_position_path    path to folder with selected measurement
          - startdate, enddate
        - part_init_config.yaml file:
          - num_part: 40000 
          - not necessarily necessary: (released species & mass, number of particles per release, layers per release, height range of released particles)
  - prepares necessary files for total column release, one FLEXPART run per day, based on specifications in options_config.yaml
    - one can list individual release positions/times in config file or read from files in satellites_position_path (created with create_GOSATpositions.py) 
  - creates:
    - subdirectory for each FLEXPART run release day: outdir/YYYY_MM/Release_YYYYMMDD/
      containing part_ic.nc file specifying initital conditions for releases, based on part_init_config.yaml
    - FLEXPART options directories: outdir/YYYY_MM/config/options_YYYYMMDD/
      start and stop in COMMAND file is adapted, rest is copied form options_dummy_directory
    - specified number of pathnames directories per month, outdir/YYYY_MM/config/pathnames_i
      one pathnames file per FLEXPART run, can start runs individually (using slurm_flexpart_v11.sh) or all within one directory (using slurm_start_multiiple.sh), or mltiple directories (using Bulkstart_multiple.sh)
    - creates part_init.nc file for total column release to use as userdefined initial conditions for FLEXPART_v11 (inside config/options_YYYYMMDD/)
  - **total column release defined in create_part_init.py**


- prepare_ISruns.py
  - run with: *python prepare_ISruns.py --config_path configs/options_config_IS.yaml*
  - **adapt before running**:
      - adapt 
        - config file: configs/options_config_RemoTeC240.yaml
          - options_dummy_path   path to '.../PK_FarewellPackage/software_levante/flexpart_v11_start_multiple/configs/options_dummy'
          - output_path 
          - input_paths         path to ERA5 data
          - available_paths     path to AVALILABLE file for ERA5 data
          - part_init_config_path
          - sim_length
          - IS_positions_dir    path to folder with selected measurement
          - startdate, enddate
        - part_init_config.yaml file:
          - num_part: 4000??
          - not necessarily necessary: (released species & mass, number of particles per release, layers per release, height range of released particles)
  - prepare necessary files for insitu release, one FLEXPART run per day, based on specifications in options_config.yaml
  - release of particles over 4h period starting at time specified in ISpositions.csv files, at given location at intake_height[magl]
  - creates same file structure as prepare_GOSATruns.py
  - creates part_init.nc file for release to use as userdefined initial conditions for FLEXPART_v11

- prepare_TCCONruns_v2.py
  - run with: *python prepare_GOSATruns.py --config_path configs/options_config_TCCON.yaml*
  - **adapt before running**:
      - adapt 
        - config file: configs/options_config_RemoTeC240.yaml
          - options_dummy_path   path to '.../PK_FarewellPackage/software_levante/flexpart_v11_start_multiple/configs/options_dummy'
          - output_path 
          - input_paths         path to ERA5 data
          - available_paths     path to AVALILABLE file for ERA5 data
          - part_init_config_path
          - sim_length
          - TCCON_measurements_path    path to folder with selected measurement
          - startdate, enddate
        - part_init_config.yaml file:
          - num_part: 40000 
          - not necessarily necessary: (released species & mass, number of particles per release, layers per release, height range of released particles)
  - prepares necessary files for total column release, one FLEXPART run per day, based on specifications in options_config.yaml
  - creates same file structure as prepare_GOSATruns.py
    - creates part_init.nc file for total column release to use as userdefined initial conditions for FLEXPART_v11 (inside config/options_YYYYMMDD/)
  - **total column release defined in create_part_init_TCCON_v2.py**


## Starting the FLAXPART runs:
- slurm_flexpart_v11.sh
  - start one FLEXPART run based on pathnames file
  - **adapt before running**:
    - adapt error and ouput filepath for the slurm jobs
    - adapt conda environment and path to FLEXPART executable
    - set pathnames path
- slurm_start_multiple.sh
  - start multiple FLEXPART runs (within one slurm job), each based on pathnames file in specified directory
  - **adapt before running**:
    - adapt error and ouput filepath for the slurm jobs
    - adapt conda environment and path to FLEXPART executable
    - enable pathnames directory path, disable line for use with Bulkstart_multiple.sh  
  - creates log.txt file for each FLEAXPER run, saved in respective output directory
- Bulkstart_multiple.sh
  - start multiple slurm jobs with FLEXPART runs from multiple pathnames directories (all FLEXPART runs within a directory are run within one slurm job, see slurm_start_multiple.sh)
  - **adapt before running**:
    - adapt list of months, paths to and number of pathnames directories
    - adapt error and ouput filepath for the slurm jobs
    - adapt conda environment and path to FLEXPART executable
    - make sure pathnames def for use with Bulkstart_multiple.sh is enabled in slurm_start_multiple.sh
  - creates log.txt file for each FLEAXPER run, saved in respective output directory




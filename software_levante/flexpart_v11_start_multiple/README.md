# flexpart_v11_start_multiple

Necessary preparations for Flexpart runs:
create_GOSATpositions.py
- creates a .csv file with all GOSAT measurements inside the selected region for each day, used RemoTeC v2.4.0 for my thesis, has option to Average measurements if close together (not used for my thesis)
- creates subdirectories: outdir/YYYY_MM/outfile_YYYYMMDD.csv, files contain columns: time,latitude,longitude,xco2,xco2_err
- run with python create_GOSATpositions.py --gosat_config_path configs/gosat_config_RemoTeC240.yaml
create_ISpositions.py
- creates a .csv file for in-situ measurements for each day, based on .nc file (created with selectObsPackData.py skript)
- creates subdirectories: outdir/YYYY_MM/outfile_YYYYMMDD.csv, files contain columns: time (at start of 4h mean time period) ,file,latitude,longitude,co2_val[ppm],elevation[masl],intake_height[magl]
- run with python create_ISpositions.py --config_path configs/IS_config.yaml
prepare_GOSATruns.py
- prepare necessary files for total column release, one FLEXPART run per day, based on specifications in options_config.yaml
    - can list individual release positions/times or read from GOSATsounding_position.csv files
- creates:
    - subdirectory for each FLEXPART run release day: outdir/YYYY_MM/Release_YYYYMMDD/
      containing part_ic.nc file specifying initital conditions for releases, based on part_init_config.yaml
    - FLEXPART options directories: outdir/YYYY_MM/config/options_YYYYMMDD/
      start and stop in COMMAND file is adapted, rest is copied form options_dummy_directory
    - specified number of pathnames directories per month, outdir/YYYY_MM/config/pathnames_i
      one pathnames file per FLEXPART run, can start runs individually (using slurm_flexpart_v11.sh) or all within one directory (using slurm_start_multiiple.sh), or mltiple directories (using Bulkstart_multiple.sh)
prepare_ISruns.py
- prepare necessary files for insitu release, one FLEXPART run per day, based on specifications in options_config.yaml
- release of particles over 4h period starting at time specified in ISpositions.csv files, at given location at intake_height[magl]
- creates same file structure as prepare_GOSATruns.py
create_part_init.py
- create part_init.nc file for total column release to use as userdefined initial conditions for FLEXPART_v11
- total column release based on config/part_init_config.yaml, used in prepare_GOSATruns.py 

Starting the FLAXPART runs:
slurm_flexpart_v11.sh
- start one FLEXPART run based on pathnames file
slurm_start_multiple.sh
- start multiple FLEXPART runs (within one slurm job), each based on pathnames file in specified directory
- creates log.txt file for each FLEAXPER run, saved in respective output directory
Bulkstart_multiple.sh
- start multiple slurm jobs with FLEXPART runs from multiple pathnames directories (all FLEXPART runs within a directory are run within one slurm job, see slurm_start_multiple.sh)
- creates log.txt file for each FLEAXPER run, saved in respective output directory



# The following things have to be adapted before using the scripts:
create_GOSATpositions.py:
- adapt config file, e.g. configs/gosat_config_RemoTeC240.yaml
create_ISpositions.py
- adapt config file, e.g. configs/IS_config.yaml
prepare_GOSATruns.py
- adapt configs/options_config_RemoTeC240.yaml
- adapt configs/part_init_config.yaml (released species & mass, number of particles per release, layers per release, height range of released particles)
- adapt configs/options_dummy/OUTGRID, output grid for FLEXPART footprints
- optionally: Other options_dummy/COMMAND file options, adapt particle output variables options_dummy/PARTOPTIONS, ...
create_part_init.py
- adapt configs/part_init_config.yaml
prepare_ISruns.py
- adapt configs/options_config_IS.yaml
- adapt configs/part_init_config.yaml (released species & mass, number of particles per release, layers per release)
- adapt configs/options_dummy/OUTGRID, output grid for FLEXPART footprints
- optionally: Other options_dummy/COMMAND file options, adapt particle output variables options_dummy/PARTOPTIONS, ...
slurm_flexpart_v11.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment and path to FLEXPART executable
- set pathnames path
slurm_start_multiple.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment and path to FLEXPART executable
- enable pathnames directory path, disable line for use with Bulkstart_multiple.sh
Bulkstart_multiple.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment and path to FLEXPART executable
- make sure pathnames def for use with Bulkstart_multiple.sh is enabled in slurm_start_multiple.sh
- adapt list of months, paths to and number of pathnames directories
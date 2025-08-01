#!/bin/bash
#SBATCH --job-name=prepare_flexpart  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=60:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=slurm/prepare_flexpart.o%j    # File name for standard output
#SBATCH --error=slurm/prepare_flexpart.e%j     # File name for standard error output

# skript to convert levante ERA5 data for use with FLEXPART_v11
# adapt slurm output directory, path to flex_extract installation, path to config file, path to control file

# to start this script you need to have metview installed (conda install metview -c conda-forge, conda install metview-python -c conda-forge)
# start this script with: sbatch prepare_flexpart.sh
eval "$(conda shell.bash hook)"     # activate conda env
conda activate flex

# paths and variables
START_DATE=20110202
END_DATE=20110301
CONFIG_PATH="/work/bb1170/RUN/b382762/software/ERA5download/convertERA5/convertERA5/config/config_2011.yaml"
CONTROLFILE_PATH=/work/bb1170/RUN/b382762/data/ERA5_daten/test/2011_02/CONTROL
OUTPUT_PATH=/work/bb1170/RUN/b382762/data/ERA5_daten/test/2011_02

# echo $START_DATE

# metview conversion
echo "starting metview"
echo $CONFIG_PATH
python convert_era5_dkrz_ml_v6.py --config_path $CONFIG_PATH

# TODO add check to see if restartAnog folder is empty, otherwise need to rerun the metview skript

# flex_extract preprocessing
cd /work/bb1170/RUN/b382762/software/flex_extract/Source/Python/Mods
pwd
echo "starting prepare_flexpart.py"
python prepare_flexpart.py --start_date $START_DATE --end_date $END_DATE --controlfile $CONTROLFILE_PATH --inputdir "$OUTPUT_PATH"/temp --ppid 660332
# TODO read dates and controlfile path from config file

# resorting EA files
cd /work/bb1170/RUN/b382762/software/ERA5download/convertERA5/convertERA5/
echo "starting resorting EA files"
sbatch --export=EAfolderpath="$OUTPUT_PATH"/temp/,outfile="$OUTPUT_PATH/",year='',month='',day='' Bulkstart_prepare_EA_files.sh

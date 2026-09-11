#!/bin/bash
#SBATCH --job-name=Start_ConvertEra5  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=10:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=/work/bb1170/RUN/b382762/data/ERA5_daten/convertERA5_test/slurm/Start_ConvertEra5.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/ERA5_daten/convertERA5_test/slurm/Start_ConvertEra5.e%j     # File name for standard error output

#to start this script you need to have metview installed (conda install metview -c conda-forge, conda install metview-python -c conda-forge)
#start this script with: sbatch Start_ConvertEra5.sh
eval "$(conda shell.bash hook)"     # activate conda env
conda activate metview

python convert_era5_dkrz_ml_v6.py --config_path config/config_test.yaml

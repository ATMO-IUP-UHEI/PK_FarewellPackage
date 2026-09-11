#!/bin/bash
#SBATCH --job-name=run_inversion     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=16:00:00       # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=150G                # compute: 0 shared: 235G (for 1.5 years/ 80x80 cov matrix 150G is enough)
#SBATCH --output=slurm/run_inversion.o%j    # File name for standard output
#SBATCH --error=slurm/run_inversion.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate pyinverse

CONFIG_PATH='/work/bb1170/RUN/b383736/software/test_PK/PK_FarewellPackage/software_levante/NA_inversion/inversion_config/config_TM5_v2_2.yaml'
echo $CONFIG_PATH
python -u run_inversion_v2.py --config "$CONFIG_PATH"
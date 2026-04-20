#!/bin/bash
#SBATCH --job-name=prepare_inversion     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=60:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=200G
#SBATCH --output=slurm/prepare_inversion.o%j    # File name for standard output
#SBATCH --error=slurm/prepare_inversion.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate inversion

python prepare_inversion.py
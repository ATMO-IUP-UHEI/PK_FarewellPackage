#!/bin/bash
#SBATCH --job-name=prep_CAMS     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=30:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=235G
#SBATCH --output=slurm/prep_CAMS.o%j    # File name for standard output
#SBATCH --error=slurm/prep_CAMS.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate inversion

python prep_CAMS_bg.py 
#!/bin/bash
#SBATCH --job-name=prepare_inversion     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=10:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=40G
#SBATCH --output=/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/slurm/prepare_inversion.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/slurm/prepare_inversion.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate pyinverse

python prepare_inversion.py
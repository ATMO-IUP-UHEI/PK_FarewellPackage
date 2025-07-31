#!/bin/bash
#SBATCH --job-name=flexpart_v11     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=48:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=50G
#SBATCH --output=slurm/flexpart_v11.o%j    # File name for standard output
#SBATCH --error=slurm/flexpart_v11.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate metview

PATHNAMES_PATH="/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/insitu/2010_06/config/pathnames/pathnames_20100601"
srun /work/bb1170/RUN/b382762/software/flexpart_v11/src/FLEXPART_ETA $PATHNAMES_PATH 
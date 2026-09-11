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
conda activate inversion

PATHNAMES_PATH="/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/RemoTeCv240/2023_02/config/pathnames_0/pathnames_20230227"
echo "$PATHNAMES_PATH"
OUTPUT_PATH=$(sed -n '2p' $PATHNAMES_PATH)
echo $OUTPUT_PATH
# run flexpart 
echo "$PATHNAMES_PATH" 
echo $(sed -n '3p' $PATHNAMES_PATH)
srun  /work/bb1170/RUN/b383736/software/flexpart/src/FLEXPART_ETA $PATHNAMES_PATH > "${OUTPUT_PATH}/log.txt" 

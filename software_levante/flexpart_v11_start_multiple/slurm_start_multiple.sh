#!/bin/bash
#SBATCH --job-name=flexpart_v11_multiple     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1              # Specify number of CPUs per task
#SBATCH --time=60:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=50G
#SBATCH --output=slurm/flexpart_v11_multiple.o%j    # File name for standard output
#SBATCH --error=slurm/flexpart_v11_multiple.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate flex

# for running this skript with one pathnames directory
# PATHNAMES_DIR="/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/RemoTeCv240/2010_06/config/pathnames_0" #path to directory containing multiple pathnames files

# for running Bulkstart_multiple.submit
PATHNAMES_DIR=$1

for PATHNAMES_PATH in "$PATHNAMES_DIR"/*
do
    # get output path from pathnames file
    OUTPUT_PATH=$(sed -n '2p' $PATHNAMES_PATH)
    echo $OUTPUT_PATH
    # run flexpart 
    echo "$PATHNAMES_PATH" 
    echo $(sed -n '3p' $PATHNAMES_PATH)
    srun ./src/FLEXPART_ETA $PATHNAMES_PATH > "${OUTPUT_PATH}/log.txt"     # run the command with path to pathnames directory, save output in log.txt file in respective output directory
done

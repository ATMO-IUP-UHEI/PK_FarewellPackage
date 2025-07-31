#!/bin/bash
#SBATCH --job-name=bulkstart     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1              # Specify number of CPUs per task
#SBATCH --time=1:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=50G
#SBATCH --output=slurm/start.o%j    # File name for standard output
#SBATCH --error=slurm/start.e%j     # File name for standard error output

# run with sbatch Bulkstart_multiple.sh
# list of months, format XX
for month in 02
do
    for i in 0 1 2  # list number of pathnames directories
    do
        # gosat releases
        pathnames_dir=/work/bb1170/RUN/b382762/data/Flexpart11_invSetup_final/RemoTeCv240/2011_${month}/config/pathnames_${i}
        echo starting $pathnames_dir
        sbatch slurm_start_multiple.submit $pathnames_dir
        # insitu releases
        pathnames_dir=/work/bb1170/RUN/b382762/data/Flexpart11_invSetup_final/insitu/2011_${month}/config/pathnames_${i}
        echo starting $pathnames_dir
        sbatch slurm_start_multiple.submit $pathnames_dir
    done
    
done

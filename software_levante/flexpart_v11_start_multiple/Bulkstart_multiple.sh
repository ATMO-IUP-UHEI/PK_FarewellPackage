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
for month in 01 02 03 # 01 02 03 04 05 06 07 08 09 10 11 12
do
    for i in 0  # list number of pathnames directories
    do
        # gosat releases
        pathnames_dir=/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/RemoTeCv240/2022_${month}/config/pathnames_${i}
        echo starting $pathnames_dir
        sbatch slurm_start_multiple.sh $pathnames_dir
        # TCCON releases
        pathnames_dir=/work/bb1170/RUN/b383736/data/Flexpart_2021/Flexpart/TCCON/2022_${month}/config/pathnames
        echo starting $pathnames_dir
        sbatch slurm_start_multiple.sh $pathnames_dir
        # insitu releases
        pathnames_dir=/work/bb1170/RUN/b383736/data/PK_Flexpart/2months/insitu/2010_${month}/config/pathnames
        echo starting $pathnames_dir
        sbatch slurm_start_multiple.sh $pathnames_dir
    done
    
done

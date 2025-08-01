#!/bin/bash
#SBATCH --job-name=merging  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=03:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=/work/bb1170/RUN/b382762/data/ERA5_daten/convertERA5_test/slurm/merge/merge.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/ERA5_daten/convertERA5_test/slurm/merge/merge.e%j     # File name for standard error output


ulimit -S -s 102400

export CDO_FILE_SUFFIX="NULL"

echo Merging data into $outfilemerge
echo $filelist $outfilemerge
grib_copy -B'date:i asc,time:i asc,level:i asc' $filelist $outfilemerge

for singlefile in $filelistrm
do
    rm $singlefile
    # echo removing $singlefile
done

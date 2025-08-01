#!/bin/bash
#SBATCH --job-name=merging  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=03:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/slurm_erroroutput_v6/merge.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/slurm_erroroutput_v6/merge.e%j     # File name for standard error output


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

# not used anymore, needed for grib2 convertion
#if [[ $filelist == *"FCOG_acc_"* ]]
#then
#    echo deleting extra files
#    for singlefile in $filelist
#    do
#        singlefile2=$( echo $singlefile | cut -f1 -d'.')
#        rm ${singlefile2}.grib2
#        rm ${singlefile2}0.grib2
#    done
#fi
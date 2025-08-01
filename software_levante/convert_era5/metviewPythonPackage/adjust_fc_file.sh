#!/bin/bash
#SBATCH --job-name=$processID  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=03:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/2010_01_v6_shifted_FCOGtest/new/slurm/adjust_FCOG.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/2010_01_v6_shifted_FCOGtest/new/slurm/adjust_FCOG.e%j     # File name for standard error output

# run this script for a single file with sbatch adjust_fc_file.sh and paths:
#-fc_file: absolute path to forecast file: 
#-outfile: absolute path to output file: 
#-processID: needed for convert_era5_dkrz_ml_v6.py skript
# example:
# fc_file=/my/path/to/fcfilefolder/FCOG_acc_SLE5sf12_1H_2010-01-01_176.grb
# outfile=/my/path/to/outputfolder/FCOG_acc_SLE5sf12_1H_2010-01-01_176.grb

# or within a python script with 
# subprocess.run(["sbatch",f"--export=fc_file={outfile},outfile={outfile}",f"--job-name={processID}","adjust_fc_file.sh"])


# echo Adapt order of messages in $fc_file and save it into $outfile
grib_set -w step=13 -s step=1,stepRange=0-1,time=1800 $fc_file ${outfile}_tempfile1
grib_set -w step=14 -s step=2,stepRange=1-2,time=1800 ${outfile}_tempfile1 ${outfile}_tempfile2
grib_set -w step=15 -s step=3,stepRange=2-3,time=1800 ${outfile}_tempfile2 ${outfile}_tempfile3
grib_set -w step=16 -s step=4,stepRange=3-4,time=1800 ${outfile}_tempfile3 ${outfile}_tempfile4
grib_set -w step=17 -s step=5,stepRange=4-5,time=1800 ${outfile}_tempfile4 ${outfile}_tempfile5
grib_set -w step=18 -s step=6,stepRange=5-6,time=1800 ${outfile}_tempfile5 ${outfile}_tempfile6
grib_set -w step=19 -s step=7,stepRange=6-7,time=1800 ${outfile}_tempfile6 ${outfile}_tempfile7
grib_set -w step=20 -s step=8,stepRange=7-8,time=1800 ${outfile}_tempfile7 ${outfile}_tempfile8
grib_set -w step=21 -s step=9,stepRange=8-9,time=1800 ${outfile}_tempfile8 ${outfile}_tempfile9
grib_set -w step=22 -s step=10,stepRange=9-10,time=1800 ${outfile}_tempfile9 ${outfile}_tempfile10
grib_set -w step=23 -s step=11,stepRange=10-11,time=1800 ${outfile}_tempfile10 ${outfile}_tempfile11
grib_set -w step=24 -s step=12,stepRange=11-12,time=1800 ${outfile}_tempfile11 $outfile

#delete temp files
for num in 1 2 3 4 5 6 7 8 9 10 11
do
    rm ${outfile}_tempfile$num
done
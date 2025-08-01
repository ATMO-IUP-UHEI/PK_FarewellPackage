#!/bin/bash
#SBATCH --job-name=ReorderEAfile  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=30:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=slurm/bulk_resort.o%j    # File name for standard output
#SBATCH --error=slurm/bulk_resort.e%j     # File name for standard error output


# run this script for a single file do
# sbatch --export=EAfolderpath=PathToFolderwithEAfiles,outfile=PathToOutputFolder,year=YearofEAfiles(YY),month=MonthOfEAfiles(MM),day=DayOfEAfiles(DD) Bulkstart_prepare_EA_files.sh
# to run the script for all days within (a month and) a year, leave the day (and the month) empty e.g.
# sbatch --export=EAfolderpath=PathToFolderwithEAfiles,outfile=PathToOutputFolder,year=YearofEAfiles,month=MonthOfEAfiles,day='' Bulkstart_prepare_EA_files.sh
# sbatch --export=EAfolderpath=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/2010/2010_06_preprocessing,outfile=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/2010/2010_06,year=2010,month=6,day='' Bulkstart_prepare_EA_files.sh

# to run this script for multiple days, month... 
# for day in 1 2 3 4 5 6; do sbatch --export=EAfolderpath=PathToFolderwithEAfiles,outfile=PathToOutputFolder,year=YearofEAfiles,month=MonthOfEAfiles,day=$day BulkStart_prepare_EA_files.sh; done

#THEREBY:
#-PathToOutputFolder: absolute path to folder where the output file should be stored in (differnt from PathToEAfileFolder!): /my/path/to/folder/
#-PathToFolderwithEAfiles:  absolute path to folder with the EAfiles: /my/path/to/EAfilefolder/
#-year: last two digits of year e.g. 10 for 2010
#-month: month as MM e.g. 05
#-day: day as DD e.g. 01

ulimit -S -s 102400

export CDO_FILE_SUFFIX="NULL"

echo Reorder all EA${year}${month}${day}*

echo ${EAfolderpath}/EA${year}${month}${day}*
EAfile_list=${EAfolderpath}/EA${year}${month}${day}*

for EAfile in $EAfile_list
do
	echo reordering $EAfile

	prepare_EA_files.sh $EAfile ${outfile}
	wait
done

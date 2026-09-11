#!/bin/bash
#SBATCH --job-name=ReorderEAfile  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=03:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=slurm/resort/resort.o%j    # File name for standard output
#SBATCH --error=slurm/resort/resort.e%j     # File name for standard error output


# run this script for a single file do
# sbatch prepare_EA_files.sh PathToEAfile PathToOutputFolder

# to run this script for all EA* files in a subdirectory and save output in a new directory 
# (please replace <PathToEAfileFolder> and <PathToOutputFolder>)
# for file in <PathToEAfileFolder>/EA*; do sbatch prepare_EA_files.sh $file <PathToOutputFolder>; done

# or within a python script with (not tested yet)
# subprocess.run(["sbatch","prepare_EA_files.sh","PathToEAfile","PathToOutputFolder"])

#THEREBY:
#-PathToEAfile: absolute path to EA file, without .grb ending: /my/path/to/EAfilefolder/EA10011500
#-PathToOutputFolder: absolute path to folder where the output file should be stored in (differnt from PathToEAfileFolder!): /my/path/to/folder/
#-PathToEAfileFolder:  absolute path to folder with the EAfiles: /my/path/to/EAfilefolder

ulimit -S -s 102400

export CDO_FILE_SUFFIX="NULL"

EAfile=$1
outfile=$2

echo Adapt order of messages in $EAfile and save it into $outfile

EAfileName="${EAfile:0-10}" #file name of the EA file is assumed to have 10 digets (= last 10 digets of provided path)

#select the individual messages in the correct order and save individually in temp files
grib_copy -w levtype=ml $EAfile ${outfile}${EAfileName}tempfile1
grib_copy -w param=134 $EAfile ${outfile}${EAfileName}tempfile2
grib_copy -w param=142 $EAfile ${outfile}${EAfileName}tempfile3
grib_copy -w param=143 $EAfile ${outfile}${EAfileName}tempfile4
grib_copy -w param=146 $EAfile ${outfile}${EAfileName}tempfile5
grib_copy -w param=180 $EAfile ${outfile}${EAfileName}tempfile6
grib_copy -w param=181 $EAfile ${outfile}${EAfileName}tempfile7
grib_copy -w param=176 $EAfile ${outfile}${EAfileName}tempfile8
grib_copy -w param=141 $EAfile ${outfile}${EAfileName}tempfile9
grib_copy -w param=151 $EAfile ${outfile}${EAfileName}tempfile10
grib_copy -w param=164 $EAfile ${outfile}${EAfileName}tempfile11
grib_copy -w param=165 $EAfile ${outfile}${EAfileName}tempfile12
grib_copy -w param=166 $EAfile ${outfile}${EAfileName}tempfile13
grib_copy -w param=228004 $EAfile ${outfile}${EAfileName}tempfile14
grib_copy -w param=168 $EAfile ${outfile}${EAfileName}tempfile15
grib_copy -w param=129 $EAfile ${outfile}${EAfileName}tempfile16
grib_copy -w param=172 $EAfile ${outfile}${EAfileName}tempfile17
grib_copy -w param=160 $EAfile ${outfile}${EAfileName}tempfile18
grib_copy -w param=27 $EAfile ${outfile}${EAfileName}tempfile19
grib_copy -w param=28 $EAfile ${outfile}${EAfileName}tempfile20
grib_copy -w param=244 $EAfile ${outfile}${EAfileName}tempfile21
#merge all temp files in correct order
grib_copy ${outfile}${EAfileName}tempfile1 ${outfile}${EAfileName}tempfile2 ${outfile}${EAfileName}tempfile3 ${outfile}${EAfileName}tempfile4 ${outfile}${EAfileName}tempfile5 ${outfile}${EAfileName}tempfile6 ${outfile}${EAfileName}tempfile7 ${outfile}${EAfileName}tempfile8 ${outfile}${EAfileName}tempfile9 ${outfile}${EAfileName}tempfile10 ${outfile}${EAfileName}tempfile11 ${outfile}${EAfileName}tempfile12 ${outfile}${EAfileName}tempfile13 ${outfile}${EAfileName}tempfile14 ${outfile}${EAfileName}tempfile15 ${outfile}${EAfileName}tempfile16 ${outfile}${EAfileName}tempfile17 ${outfile}${EAfileName}tempfile18 ${outfile}${EAfileName}tempfile19 ${outfile}${EAfileName}tempfile20 ${outfile}${EAfileName}tempfile21 $outfile${EAfileName}

#delete temp files
for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
do
    rm ${outfile}${EAfileName}tempfile$num
done
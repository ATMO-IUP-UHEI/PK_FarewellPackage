Created by Eva-Marie Metz 23.04.2024, adapted by Pernilla Kühn

Program to convert ERA5 spherical data on levante to files usable by Flexextract and FLEXPART. It replaces the download of ERA5 data performed with flex_extract via the command get_mars_data. THe output data needs to be preprocessed with the prepare_flexpart script of felx_extract.

The following scripts are part of the program:
convert_era5_dkrz_ml_v6.py
- Main script to convert data for individual parameters or file format needed by flex_extract (combine the daily data for necessary parameters), based on config file
Start_ConvertEra5.sh
- start main script with sbatch Start_ConvertEra5.sh
prepare_flexpart.sh
- start all steps to prepare data from FLEXPART run with sbatch prepare_flexpart.sh, need flex_extract installation
    (includes Start_ConvertEra5.sh, start_flex_extract_prepareflexpart.sh, Bulkstart_prepare_EA_files.sh)
sh_merge.sh
- merges the individual files for each day and parameter together in 3-day files with all parameters of one type 
convert_ANOG__ML.py
- the data for the ANOG__ML file type of flex_extract is converted parallel in the background jobs. This script is used for the background jobs
ConvertANOG__ML.sh
- start the background jobs
adjust_fc_file.sh
- adjust times in forecast fields for first day of the month (due to inconsistencies in file format)
check_status.sh
- check that all background jobs are completed before merging ANOG__ML or FCOG_acc_SL files
prepare_EA_files.sh
- adapt order of parameters in EA files produced with start_flex_extract_prepareflexpart.sh for use with FLEXPART
- for single file or all files in a subdirectory
Bulkstart_prepare_EA_files.sh
- start prepare_EA_files.sh skript for all files in a day/month/year simultaneously

The following things have to be adapted before using the scripts:
config/config_test.yaml
- set time period, output path, resolution, area, control_outpath if use for flepxart
- set individual parameters/ flexpart files that should be downloaded
- for use with prepare_flexpart.sh, adapt CONTROL dummy 
Start_ConvertEra5.sh:
- adapt error and ouput filepath for the slurm jobs
- adapt config filepath and config file content
- adapt conda environment
prepare_flexpart.sh:
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment
- adapt listed paths and variables
- adapt flex_extract source path
- adapt convertERA5 source path
sh_merge.sh:
- adapt error and ouput filepath for the slurm jobs
ConvertANOG__ML.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment
adjust_fc_file.sh
- adapt error and ouput filepath for the slurm jobs
check_status.sh!!!
- adapt userID
prepare_EA_files.sh:
- adapt error and ouput filepath for the slurm jobs
Bulkstart_prepare_EA_files.sh
- adapt error and ouput filepath for the slurm jobs

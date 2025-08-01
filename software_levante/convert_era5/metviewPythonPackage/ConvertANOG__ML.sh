#!/bin/bash
#SBATCH --job-name=$processID  # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --nodes=1              # Specify number of nodes
#SBATCH --mem=24000            # Specify memory to be used for job (MB)
#SBATCH --time=03:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --output=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/slurm_erroroutput_v6/ANOG.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/ERA5_daten/northAmerica/slurm_erroroutput_v6/ANOG.e%j     # File name for standard error output

#this script is part of convert_era5_dkrz_ml_v5.py
#conda deactivate
source /sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/etc/profile.d/conda.sh
conda deactivate
conda activate metview
  
echo $CONDA_PREFIX

python convert_ANOG__ML.py $filepath $inter_res $outfile $par $resol $borders0 $borders1 $borders2 $borders3

# delete restart file
rm $restartfile
echo deleted restart file $restartfile

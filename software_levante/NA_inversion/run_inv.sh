#!/bin/bash
#SBATCH --job-name=run_inversion     # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1     # Specify number of CPUs per task
#SBATCH --time=30:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=bb1170       # Charge resources on this project account
#SBATCH --mem=235G
#SBATCH --output=/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/slurm/run_inversion.o%j    # File name for standard output
#SBATCH --error=/work/bb1170/RUN/b382762/data/FarewellPackage_test/Flexpart/slurm/run_inversion.e%j     # File name for standard error output

eval "$(conda shell.bash hook)"     # activate conda env
conda activate pyinverse

CONFIG_PATH='/work/bb1170/RUN/b382762/master_farewell_package/software_levante/NA_inversion/inversion_config/config_test.yaml'
echo $CONFIG_PATH
python run_inversion.py --config "$CONFIG_PATH"
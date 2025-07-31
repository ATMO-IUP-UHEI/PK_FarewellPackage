# inversion of NA CO2 fluxes using FLEXPART

prepare_inversion.py
- necessary preprocessing 
    - calculate background values for each measurement, 
    - prepare TM5-4DVar prior and RemoTeC flux to use as reference, 
    - prepare footprints (add diurnal cycle scaling, sum to weekly 2x2 resolution) 
    - optionally: calculate fraction of remaining particles at simulation end
- see function definitions for detail
- should be started using prep_inv.sh for large number of measurements
prep_inv.sh
- script to start prepare_inversion.py

run_inversion.py
- main skript to run the inversions, based on configurations given in inversion_config/config.yaml
run_inv.sh
- script to start run_inversion.py

analyze_inversion.py
- calculate modeled measurement values 
    - using our posterior fluxes from the performed inverisons and weekly footprints
    - from high resolution footprints (1x1, hourly) and TM5-4DVar fluxes
- get mean and range of measurements/modeled values for each week
- plot 
    - timeseries of prior, posterior and TM5-4DVar reference for 2x2 inner domain (or 4x4 inversion entire domain)
    - difference histogramm between modeled measurement values and measurements
- see function definitions for detail
analyze_inv.sh
- script to start analyze_inversion.py

...
TODO


The following things have to be adapted before using the scripts:
prepare_inversion.py
- adapt paths in main function
prep_inv.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment
run_inversion.py
- adapt config file
run_inv.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment and config path
analyze_inversion.py
- select which functions should be run, adapt necessary paths, list of measurement error values, footprint variable name etc
- 
- 
analyze_inv.sh
- adapt error and ouput filepath for the slurm jobs
- adapt conda environment
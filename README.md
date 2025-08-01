# master_farewell_package

Folder structure
- Figures – all final figures included in my master’s thesis
- Software_levante
    - convertERA5 packaage to convert ERA5 data to lat/lon coordinates, can be used for individual parameters or to get all data for a FLEXPART run (still needs to be preprocessed using starting prepare_flexpart.py from flex_extract, and resort produced EA files)
    - flexpart_v11_start_multiple: skripts to prepare and start FLEXPART runs for GOSAT retrievals and In-situ measurements
    - NA_inversion: Skripts to prepare, run and first analysis of CO2 flux inversion
- Software_atmo
    Prepare_TM5-4DVar data diurnal cycle scaling factors, molefraction fields and preprocess ObsPack data (see README for details)

see folder/README for more details

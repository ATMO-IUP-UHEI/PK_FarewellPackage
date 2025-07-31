# Farewell package Pernilla Kühn

prep_TM5-4DVar.py
- Caclulate diurnal cycle scaling factors based on 3hrly TM5-4DVar molefraction fields
- Cut desired region, calculate pressure at TM5-4DVar layer boundaries, calculate xco2 from TM5-4DVar molefraction fields
selectOnsPackData.py
- preprocessing of obspack_co2_1_GLOBALVIEWplus_v5.0_2019-08-12 dataset, 
    Selects surface & tower measurements within defined region & time period, only uses highest tower inlet, CarbonTracker assimilation flag CT_assim=1, caclulates 4h mean


The following things have to be adapted before using the scripts:
prep_TM5-4DVar.py
- adapt desired region, time period, output paths
- select desired datasets for molefraction fields
selectOnsPackData.py
- adapt desired region, time period, output paths

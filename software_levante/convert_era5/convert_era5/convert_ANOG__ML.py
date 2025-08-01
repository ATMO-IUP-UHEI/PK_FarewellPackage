# Script to automate conversion of ERA5 (spectral) files on DKRZ Levante to flexextract download files
import datetime as dt
import subprocess
import argparse
import glob
import os
import pandas as pd
import time
import metview as mv

           
#this Script is not supposed to be run like this but started by ConvertANOG__ML.sh which is part of convert_era5_dkrz_ml_v5.py

def parse_arguments_ANOG():
    parser = argparse.ArgumentParser(description="Script to automate conversion of ERA5 spectral files into flex_extract download files (downloaded with get_mars_data e.g.'python get_mars_data.py --start_date 20100110 --end_date 20100115 --controlfile CONTROL_EA5_global_coarse --inputdir 'apth/to/output/file/'). ######################################################## run this script e.g. with: python convert_era5_dkrz_ml.py '20091201' '20091211' '/work/bb1170/RUN/b381317/data/FLEXPART_ECMWF/ERA5/Levante/EA_box_20091201-20091211-Levante/' 10 -170 -90 180 90 --fe_files 'FCOG_acc_SL'")
    parser.add_argument("filepath",    help="Grib file path")
    parser.add_argument("inter_res",  help="inter res for interpolation")
    parser.add_argument("outfile",    help="Grib output file path")
    parser.add_argument("parID",    help="parameter number")
    parser.add_argument("res",    help="resolution for interpolation")
    parser.add_argument("borders", type=float, nargs='+',   help="borders as list: lat_min, long_min, lat_max, long_max")
    args = parser.parse_args()
        
    filepath = args.filepath
    res = float(args.res)
    inter_res = int(args.inter_res)
    borders = args.borders
    outfile = args.outfile
    parID = str(args.parID)
    
    return filepath, res, inter_res, borders, outfile, parID

def convert_anog_ml(filepath, res, inter_res, borders, outfile, parID):
    datagrb = mv.read(filepath)
    
    datagrb_interp = mv.read(data=datagrb,    
                            grid=[str(res),str(res)],
                            resol=inter_res,
                            accuracy=24,
                            interpolation='linear',
                            area=borders)

    mv.write(outfile + '_' + str(parID)+'.grb', datagrb_interp)


def main():
    filepath, res, inter_res, borders, outfile, parID = parse_arguments_ANOG()
    print('create ' + outfile + '_' + str(parID)+'.grb' + ' in background')
    convert_anog_ml(filepath, res, inter_res, borders, outfile, parID)
    



##########################################################################################################
if __name__ == "__main__":
    main()
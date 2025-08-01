# Script to automate conversion of ERA5 (spectral) files on DKRZ Levante to flexextract download files
import datetime as dt
import subprocess
import argparse
import glob
import os
import pandas as pd
import numpy as np
import time
import metview as mv
import click
import yaml
import ast
import itertools
import sys
from dateutil.relativedelta import relativedelta

###################################################
# run this script e.g. with 
# python convert_era5_dkrz_ml_v5.py --config_path config/config_convert_era5.yaml
###################################################


# flex_extract files, corresponding levante folders and parameters
@click.command()
@click.option("--config_path", required=True, help="path to config.yaml file")
def main(config_path):
    ymdi, ymdf, outpath, res, borders,parameter_list, dlevel,dtype,timeres,flex_extract_download, fe_files, idnum = read_config(config_path)
    # print('started main, read config')
    if flex_extract_download:
        # read paths and params from csv file
        dfFiles= pd.read_csv('config/flexextract_files.csv', sep='\t')
        # Convert the 'paramlist' column back to a list
        dfFiles['paramlist'] = dfFiles['paramlist'].apply(ast.literal_eval)
        
        # Evas skript
        processID = int(time.time()) # needed for ANOG__ML conversion
        print('process ID = '+str(processID))

        if fe_files is None:
            print('no files given, processing all files')
            fe_files = ['ANOG__ML', 'ANOG__SL', 'ANSH__SL','FCOG_acc_SL','OG_OROLSM__SL']
        else:#check that all given file names are valid
            all_valid = True
            for fe_file_i in fe_files:
                all_valid = all_valid * (fe_file_i in dfFiles.fefile.values)
            if all_valid == 0:
                print('(one of) the given flexextract file(s) name does not exist, please try again e.g. with "FCOG_acc_SL"')
                quit()
        if idnum is None:
            idnum = '000'
        
        for fe_file in fe_files:
            print('Converting '+fe_file)
            #convert individual ERA5 files from levante from spectral_complex/reduced_gg to regular_ll for all days individually
            convertData(ymdi, ymdf, outpath, res, borders, fe_file,idnum,processID,dfFiles)
            #merge data to combined files of 3-days each
            if 'ANOG__ML' in fe_file:
                Wait_for_Jobs(processID)
            print('Merging '+fe_file)
            Merge3Days(ymdi, ymdf, outpath, fe_file,dfFiles)
    else:
        print('reading config/param_files.csv')
        # read paths and params from csv file
        dfFiles_param= pd.read_csv('config/param_files.csv', sep='\t')
        # Convert the 'paramlist' column back to a list
        dfFiles_param['paramlist'] = dfFiles_param['paramlist'].apply(ast.literal_eval)
        for param in parameter_list:
            convertParam(ymdi, ymdf, outpath, res, borders, param, dlevel,dtype,timeres,dfFiles_param)
        

def convertData(ymdi, ymdf, outpath, res, borders, fe_file,idnum, processID,dfFiles):
    '''
    # convert individual ERA5 files from levante from spectral_complex or reduced_gg to regular_ll for all days individually
    # arguments:
    #            ymdi: Start date as datetime date object
    #            ymdf: End date as datetime date object
    #            fe_file: flexextract file which is aimed to be reproduced
    # returns:
    #            nothing, creates transformed daily ERA5 files for each parameter
    '''
    data_name = dfFiles[(dfFiles.fefile == fe_file)].levanteData.values[0].split('_')
    file_name = dfFiles[(dfFiles.fefile == fe_file)].levanteFile.values[0]
    
    timerange = range((ymdf-ymdi).days+1)
    if 'FCOG' in fe_file: # FCOG files start one day before and extend to one day after the time period
        timerange = range(-1,(ymdf-ymdi).days+2)
    
    for day in timerange:
        #get date of day as datetime and string
        ymdn=ymdi+dt.timedelta(days=day)
        ymdn_str = str(ymdn.year)[0:4]+str(ymdn.month).zfill(2)+str(ymdn.day).zfill(2)
        
        #filepath for normal files
        inpath=f"/pool/data/ERA5/E5/{data_name[0]}/{data_name[1]}/{data_name[2]}/yyy/" #yyy placeholfer for parameter number later
        filename=f"E5{file_name}_{ymdn}" 
        infile=inpath+filename
        
        #filepath for invariant files see /pool/data/ERA5/READMEs7README_ERA5_IV_files.txt
        inpathIV=f"/pool/data/ERA5/E5/{data_name[0]}/{data_name[1]}/IV/yyy/" #yyy placeholfer for parameter number later
        filenameIV=f"E5{file_name.replace('1H','IV')}_2000-01-01"
        infileIV=inpathIV+filenameIV
        
        #outfile
        outfile=outpath+fe_file+filename

        #convert data depending on the individual file kind
        #ansh remains spherical and only the level and time needs to be selected    
        for par in dfFiles[(dfFiles.fefile == fe_file)].paramlist.values[0]:
            print('Converting Parameter '+ str(par) + ' on ' + ymdn_str)
            if par in ['027','028','160','129','172']: #these are the invariant files
                filepath = infileIV.replace('yyy',str(par))+'_'+str(par)+'.grb'
            else:
                filepath = infile.replace('yyy',str(par))+'_'+str(par)+'.grb'

            if res > 2:#this is our choice for the intermediate grid
                inter_res = 159
            else:
                inter_res = 799

            if 'ANOG__ML' in fe_file:#start sbatch jobs for each day and paramter to run parallel in background
                print('converting paramter '+str(par)+' for the day '+ ymdn_str +' in background')
                subprocess.run(["sbatch",f"--export=filepath={filepath},resol={res},inter_res={inter_res},outfile={outfile},par={par},borders0={borders[0]},borders1={borders[1]},borders2={borders[2]},borders3={borders[3]}",f"--job-name={processID}","ConvertANOG__ML.sh"])
                
            else:
                datagrb = mv.read(filepath)

                if 'ANSH' in fe_file:#data remains in spheric coordinates
                    datagrb_interp = mv.read(data=datagrb)
                    datagrb_interp = datagrb_interp.select(level=1)
                    datagrb_interp = mv.read(data=datagrb_interp, accuracy = 24)
                else:
                    datagrb_interp = mv.read(data=datagrb,    
                        grid=[str(res),str(res)],
                        resol=inter_res,
                        accuracy=24,
                        interpolation='linear',
                        area=borders) # S,W,N,E
                #certain specifications
                #invaiant files
                if par in ['027','028','160','129','172']: #change date from 2000-01-01 to correct date

                    datagrb_interp = mv.grib_set(datagrb_interp,['date',int(ymdn_str),'time',0000])

                    if par == '172' or par == '129':#this parameter is needed at every time but not provided like that
                        datagrb_interp0 = datagrb_interp
                        for i in range(100, 2400, 100):
                            datagrb_interp = mv.merge(datagrb_interp,mv.grib_set(datagrb_interp0,['time',i]))
                
                if 'OROLSM' in fe_file:#only time 0000 needed
                    datagrb_interp = datagrb_interp.select(time=0000)
                
                mv.write(outfile + '_' + str(par)+'.grb', datagrb_interp)

def convertParam(ymdi, ymdf, outpath, res, borders, paramid, dlevel,dtype,timeres,dfFiles_param): #TODO processid und idnum?
    '''
    # convert individual parameter from levante from spectral_complex or reduced_gg to regular_ll for all days individually
    # arguments:
    #            ymdi: Start date as datetime date object
    #            ymdf: End date as datetime date object
    #            res: 
    #            paramid: parameter that is to be converted
    # returns:
    #            nothing, creates transformed daily ERA5 files for each parameter
    '''
    # check if dtype and dlevel given in config fle or not
    if dlevel == [None]:
        dlevel=['sf','ml','pl']
    if dtype == [None]:
        dtype=['an','fc']
    # list of folder strings 
    # print(dtype)
    if timeres=='IV':
        folder_str = [f'{dl}/{dt}/IV' for dl, dt in list(itertools.product(dlevel, dtype))]
    else:
        folder_str = [f'{dl}/{dt}/timeres' for dl, dt in list(itertools.product(dlevel, dtype))]
    # print(folder_str)
    # get indexes of releant row of the dfFiles_param dataframe
    i_temp = dfFiles_param.index.where(dfFiles_param.levanteFolders.str.startswith(tuple(folder_str))).dropna().astype(int).to_list()
    # print(i_temp)
    if len(i_temp)>1:
        i_new=[]
        for i in range(0,len(i_temp)):
            if paramid in dfFiles_param.paramlist[i_temp[i]]: 
                i_new.append(i_temp[i])
        i_temp=i_new
        # print(i_temp)
        if len(i_temp)>1:
            temp_str=[dfFiles_param.levanteFolders[i] for i in i_temp]
            print(f'ERROR: parameter {paramid} found more than once: {temp_str}, please specify desired data type and level in config file')
            sys.exit(1)
    if len(i_temp)==0:
        print(f'ERROR: parameter {paramid} not availible for chosen data type, level and time resolution')
        sys.exit(1)
    if not paramid in dfFiles_param.paramlist[i_temp[0]]:
        print(f'ERROR: parameter {paramid} not availible for chosen data type, level and time resolution')
        sys.exit(1)

    data_name=dfFiles_param.levanteFolders[i_temp[0]].split('/')
    data_name[-1]=timeres                                          #set time res
    file_name=dfFiles_param.levanteFile[i_temp[0]][0:5]+timeres    #set time res
    # print(data_name,file_name)
    
    if timeres=='1H':
        timerange = range((ymdf-ymdi).days+1)
    if timeres=='1D':
        timerange = range((ymdf.year-ymdi.year)*12+(ymdf.month-ymdi.month)+1)
    if timeres=='1M':
        timerange = range((ymdf.year-ymdi.year)+1)

    for timestep in timerange:
        if timeres=='1H':
            #get date of day as datetime and string
            ymdn=ymdi+relativedelta(days=timestep)
            ymdn_str = ymdn.strftime("%Y-%m-%d")
        if timeres=='1D':
            ymdn=ymdi+relativedelta(months=timestep)
            ymdn_str = ymdn.strftime("%Y-%m")
        if timeres=='1M':
            ymdn=ymdi+relativedelta(years=timestep)
            ymdn_str = ymdn.strftime("%Y")
        
        #filepath for normal files
        param_str = f"{paramid:03d}"
        inpath=f"/pool/data/ERA5/E5/{data_name[0]}/{data_name[1]}/{data_name[2]}/{param_str}/" 
        
        if data_name[-1]=='IV': # for invariant files
            filename=f"E5{file_name}_2000-01-01_{param_str}.grb"
        else:
            filename=f"E5{file_name}_{ymdn_str}_{param_str}.grb" 
        infile=inpath+filename
        
        #outfile
        outfile=outpath+filename
        # for invariant files, break loop if file already exists
        if data_name[-1]=='IV' and os.path.exists(outfile):
            break
        #convert data depending on the individual file kind
        print('Converting Parameter '+ param_str + ' on ' + ymdn_str)

        inter_res=np.floor((180/res*2-1)/2)
        
        print(f'inter_res: {inter_res}')
        datagrb = mv.read(infile)
        datagrb_interp = mv.read(data=datagrb,    
            grid=[str(res),str(res)],
            resol=inter_res,
            accuracy=24,
            interpolation='linear',
            area=borders)# S,W,N,E
        # invaiant files
        # TODO just once, not for every timestep for invariant files
        if data_name[2]=='IV':
            # TODO
            #change date from 2000-01-01 to correct date
            datagrb_interp = mv.grib_set(datagrb_interp,['date',int(ymdn_str),'time',0000])
        mv.write(outfile, datagrb_interp)

def Merge3Days(ymdi, ymdf, outpath, fe_file,dfFiles):
    '''
    # merge data to combined files of 3-days each
    # arguments:
    #            ymdi: Start date as datetime date object
    #            ymdf: End date as datetime date object
    #            fe_file: flexextract file which is aimed to be reproduced
    # returns:
    #            nothing, creates flex_extract ERA 5 files and deletes files created by convert data
    '''
    data_name = dfFiles[(dfFiles.fefile == fe_file)].levanteFile.values[0]
    print('Merging data')
    
    timerange = range(0,(ymdf-ymdi).days+1,3)
    date_max = ymdf
    
    if 'FCOG' in fe_file: # FCOG files start one day before and extend to one day after the time period
        timerange = range(-1,(ymdf-ymdi).days+2,3)
        date_max = ymdf+dt.timedelta(days=1)
    for day in timerange:
        ymdn3 = ymdi+dt.timedelta(days=day)
        filelist = '' #list of daily files merged into the 3-day files, for OG_OROLSM files, only the first of each three days is needed
        filelistrm = '' #list of all daily files, they will be removed
        outfilemerge = outpath+f"{fe_file}.{ymdn3.year}"+str(ymdn3.month).zfill(2)+str(ymdn3.day).zfill(2)+".557856.660332.grb"
        for param in dfFiles[(dfFiles.fefile == fe_file)].paramlist.values[0]:
            for di in range(0,3):
                ymdn3i = ymdn3+dt.timedelta(days=di)
                if ymdn3i <= date_max: 
                    filelistrm = filelistrm + f" {outpath}{fe_file}E5{data_name}_{ymdn3i}_{param}.grb"
                    #for OG_OROLSM files, only the first of each three days is needed
                    if not ('OG_OROLSM' in fe_file and di > 0):
                        filelist = filelist + f" {outpath}{fe_file}E5{data_name}_{ymdn3i}_{param}.grb"
            
        subprocess.run(["sbatch",f"--export=filelist={filelist},filelistrm={filelistrm},outfilemerge={outfilemerge}","sh_merge.sh"])

def Wait_for_Jobs(processID):
    '''
    # wait until all jobs to convert ANOG__ML files are finished
    # arguments:
    #            Joblist: list of jobs which are waited for to finish
    '''
    running = True    
    while running == True:
        time.sleep(30)
        #check_status.sh echos 'running' when job with processID in its NAME is in squeue, 'Finished' otherwise
        result = subprocess.run(["check_status.sh",str(processID)], capture_output=True, text=True)
        if 'Finished' in result.stdout:
            if 'running' in result.stdout: # just a test
                print('something went wrong in Wait_for_Jobs')
            else:
                running = False
                print('all ANOG_ML files are created')
        else:
            print('waiting for ANOG_ML files to be created, status of jobs will be checked again in 30 sec.')

def read_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    
        # get start and stop date components
        yi=int(config['startdate'][0:4])
        mi=int(config['startdate'][4:6])
        di=int(config['startdate'][6:8])

        yf=int(config['stopdate'][0:4])
        mf=int(config['stopdate'][4:6])
        df=int(config['stopdate'][6:8])

        borders = config['borders']
        
        if len(borders) != 4 or borders[0] > borders[2] or borders[1] > borders[3]:
            SystemExit('enter 4 borders in the follwoing order lat_min, long_min, lat_max, long_max')
        
        
        
        return (dt.date(yi,mi,di), dt.date(yf,mf,df), config['outpath'], config['res'], config['borders'], config['param_list'], [config['data_level']],[config['data_type']],config['timeres'],config['flex_extract'], config['fe_files'], config['idnum'])

#########################################################################################

if __name__ == '__main__':
    main()

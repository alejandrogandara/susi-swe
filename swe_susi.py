#%% Importing general libraries and functions

import numpy as np
import pandas as pd
from datetime import datetime
import configparser
import argparse
import matplotlib.pylab as plt
import seaborn as sns
import sys
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline as interS
from scipy.sparse import  diags
from scipy.sparse.linalg import  spsolve
#----------------------------------------------------------------------------
# Install non-standard packages
# !pip install xlsxwriter
# !pip install xlrd==1.2.0
# !pip install pyproj
# #--------------------------
# !pip install netCDF4
# !pip install datetime
# !pip install xlrd --upgrade
# !pip install dataframe_image
# # --------------------------
from pyproj import CRS, Transformer
import os

# Importing repositories 
if not os.path.exists('susi_2022'): 
    !git clone https://github.com/annamarilauren/susi_2022
else: print("Susi repository is already in yor computer")

if not os.path.exists('susi_SMHI'): 
    !git clone https://github.com/alejandrogandara/susi_SMHI
else: print("SMHI collect repository is already in yor computer")

# Importing local libraries

import scripts.import_param as importParam # script to import from xls to dict
import susi_SMHI.susi_SMHI as smhi

newpathsusi = 'susi_2022'
if newpathsusi not in sys.path: sys.path.insert(0,newpathsusi)
newpathsusi = 'susi_SMHI/'
if newpathsusi not in sys.path: sys.path.insert(0,newpathsusi)
# -----------------------------------------------------------------------------
from susi.susi_utils import read_FMI_weather
from susi.susi_main import Susi
from inputs.susi_para import get_susi_para
from susi.susi_io import *

print (f"SUSI is ready to run!")

# gets parameters from  file
# This script provides functionality to get data from an xmlx file. You will find a template at ./inputs/param_template.xlsx
# pFile: parameter file

from pandas import read_excel
import traceback
from datetime import datetime

def getParamFile(path):
    #this variable was planned to 
    try: pFile = read_excel(path)
    except: 
        print('file does not exists')
        return
    return pFile

def getVariablesFromFile(pFile, site):
    # var [True, False]
    # file ['file path']
    vars = {}
    pFile=pFile[(pFile['commented'] != '#') & (pFile['path'].isna())]
    pFile.reset_index()
    for ind in pFile.index:
        if pFile['type'][ind] == 'date': 
            
            date = datetime.strptime(pFile[site][ind], '%Y-%m-%d %H:%M:%S')
            globals()[pFile['var'][ind]] = date
            vars[pFile['var'][ind]] = date
            #datetime.datetime(date.year,date.month,date.day)
        else: 
            globals()[pFile['var'][ind]] = pFile[site][ind] 
            vars[pFile['var'][ind]] = pFile[site][ind]
        print('variable updated: ' + pFile['var'][ind]+' = '+ str(pFile[site][ind]))
    return vars

# gets parameters from param file and fill dictionaries in memory
def getDicInFile_original(pFile, site):
    pFile=pFile[(pFile['commented'] != '#') & (pFile['path'].notnull())]
    pFile.reset_index()

    for uvar in pFile['var'].unique():
        sp = pFile[pFile['var'] == uvar]
        print('------------------------------')
        print('dictionary: '+ uvar + '  keys:' + str(len(sp)))
        globals()[uvar]
        for ind in sp.index:
            try: 
                path = pFile['path'][ind].replace(', ','][')
                if pFile['type'][ind] == 'str': 
                    value = "'"+ str(pFile[site][ind]) + "'"
                else: value = str(pFile[site][ind])
                #exec(pFile['var'][ind] + pFile['path'][ind].replace(',','][').replace(' ', '') + ' = ' + str(pFile[site][ind]))
                line = uvar + path + ' = ' + value
                print(line)
                exec(line)
            #except: print('ERROR: '+ uvar + path + ' = ' + value + ' NOT UPDATED' )
            except Exception: traceback.print_exc()



def getDicInFile(pFile, site):
    pFile=pFile[(pFile['commented'] != '#') & (pFile['path'].notnull())]
    pFile.reset_index()

    for uvar in pFile['var'].unique():
        sp = pFile[pFile['var'] == uvar]
        print('------------------------------')
        print('dictionary: '+ uvar + '  keys:' + str(len(sp)))
        globals()[uvar]
        for ind in sp.index:
            try: 
                path = pFile['path'][ind].replace(', ','][')
                if pFile['type'][ind] == 'str': 
                    value = "'"+ str(pFile[site][ind]) + "'"
                else: value = str(pFile[site][ind])
                #exec(pFile['var'][ind] + pFile['path'][ind].replace(',','][').replace(' ', '') + ' = ' + str(pFile[site][ind]))
                line = uvar + path + ' = ' + value
                print(line)
                exec(line)
            #except: print('ERROR: '+ uvar + path + ' = ' + value + ' NOT UPDATED' )
            except Exception: traceback.print_exc()


# Name convention
# weather file:     inputs/.../weather/  {site}_weather.csv
# motti file:       inputs/.../motti/    {site}_motti_lyr_0.xlsx
# You will define the input folder later

#%%  Parameters
folderName= importParam.mkfolder('outputs/')               # output folder
graph_folder = importParam.mkfolder(f'{folderName}/graphs/')

weatherPath = 'inputs/sweden/weather/'                          # Where the weather file is
#mottipath = 'inputs/sweden/motti/'                         # folder where motti files are
mottipath = 'inputs/sweden/heureka/'                         # folder where motti files are


paramFile = getParamFile("inputs\sweden\parameters.xlsx")  # map the param file

sitesNames = paramFile.columns[9:]  # extract the site names from param file header

# Defines motti files, please follow the name convention or add it manually
data_files = {}
for key in sitesNames:
    data_files[key] = {}
    data_files[key]['para'] = f'susi_para_{key}'
    data_files[key]['weatherFile'] = f'{key}_weather.csv'
    data_files[key]['heureka_lyr_0'] = f'{key}_heureka_input_lyr_0.xlsx'
    data_files[key]['heureka_lyr_1'] = f'{key}_heureka_input_lyr_0.xlsx'        # fake, just to fill the gap
    data_files[key]['heureka_lyr_2'] = f'{key}_heureka_input_lyr_0.xlsx'        # fake, just to fill the gap
    

#%%
##### ciclo

print(pd.DataFrame(sitesNames))
siteID = slice(7,8)  #slice, remember that python counts (#,#] so use the second position + 1
#list(data_files.keys())[1:2]
site_list = list(data_files.keys())[siteID]
print('running in:')
print(site_list)

########


#%% Getting the default parameter values
for site in site_list:
    mottifile = {'path':mottipath,
                'dominant':{1: data_files[site]['heureka_lyr_0']},
                'subdominant':{0:data_files[site]['heureka_lyr_1']},
                'under':{0:data_files[site]['heureka_lyr_2']}}

    weatherData = data_files[site]['weatherFile']
    sarkaSim = paramFile.loc[paramFile.path == "['start_date']"][site].values[0]
    start_date = paramFile.loc[paramFile.path == "['start_date']"][site].values[0]
    end_date = paramFile.loc[paramFile.path == "['end_date']"][site].values[0]
    start_yr = start_date.year 
    end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25

    sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
    n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

    #age = paramFile.loc[paramFile.path == "['age']"][site].values[0]
    dominant_age = float(paramFile.loc[paramFile['var'] == "age_dominant"][site].values[0])
    #dominant_age = 85.

    ageSim = {'dominant': dominant_age * np.ones(n),
            'subdominant': 0*np.ones(n),
            'under': 0*np.ones(n)}                                                         # age of the stand in each node

    sfc =  np.ones(n, dtype=int)*3                                                            # site fertility class      
    wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat='krycklan', 
                                                                            folderName=folderName, hdomSim=None,  
                                                                            ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath = weatherPath,
                                                                            n=n)
                                                                    
    ##%% #Get the parameters to dictionary
    getDicInFile(paramFile, site)
    spara['scenario name']

    print('1111-------------------------------------------- scenarios')
    print(spara['scenario name'])
    print('1111-------------------------------------------- scenarios')

    try:
        #Collect the weather from SMHI
        workingFolder = os.getcwd()
        wfile = smhi.getWeather(site, wpara['lon'], wpara['lat'], start_date, end_date, weatherPath, workingFolder, stations_nearby= 15)
    except: 
        print(f'error in weather data for {site}')
        continue
    
    # Feed SUSI with the weather data
    forc=read_FMI_weather(0, start_date, end_date, sourcefile=wfile)
    
    from susi.susi_io import *
    #graph_output = graph_folder + '/' + weatherData[:-12] +'/' + weatherData[:-4] + '_graph.png'
    weather_fig(forc)        #Draw the weather input figure

    print('-------------------------------------------- scenarios')
    print(spara['scenario name'])
    print('-------------------------------------------- scenarios')


    try:
        print("RUNING SUSI")
        susi = Susi()
        susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'krycklan', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)

    except: 
        print(f"\nNO SUSI FOR YOUR {site}")
        continue
print(f'\n\nEOF')
# %%

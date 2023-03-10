from pandas import read_excel
import traceback
from datetime import datetime
import os
#import operator
#import json
#from functools import reduce  # forward compatibility for Python 3

# create new folder if not exists already
def mkfolder(path):
    if not os.path.exists(path):
        os.mkdir(path)
        print(f'{path} has been created, ready to go')
    else:
        print(f'{path} already exists, ready to go')
    return path


#get variables
def getVariablesFromFile(file, site, var=False):
    vars = {}
    try: p = read_excel(file)
    except: 
        print('file does not exists')
        return
    if var:
        try: p=p[(p['commented'] != '#') & (p['path'].isna()) & (p['var'] == var)]
        except: 
            print('variable not found')
            return
    else: p=p[(p['commented'] != '#') & (p['path'].isna())]
    p.reset_index()
    for ind in p.index:
        if p['type'][ind] == 'date': 
            date = datetime.strptime(p[site][ind], '%Y-%m-%d %H:%M:%S')
            globals()[p['var'][ind]] = date
            vars[p['var'][ind]] = date
            #datetime.datetime(date.year,date.month,date.day)
        else: 
            globals()[p['var'][ind]] = p[site][ind] 
            vars[p['var'][ind]] = p[site][ind]
        print('variable updated: ' + p['var'][ind]+' = '+ str(p[site][ind]))
    return vars


def getDicInFile(file, site, var=False):
    try: p = read_excel(file)
    except: 
        print('file does not exists')
        return
    #p = p.where(notnull(p), '')
    if var:
        try: p=p[(p['commented'] != '#') & (p['path'].notnull()) & (p['var'] == var)]
        except: 
            print('variable not found')
            return
    else: p=p[(p['commented'] != '#') & (p['path'].notnull())]
    p.reset_index()

    for uvar in p['var'].unique():
        sp = p[p['var'] == uvar]
        print('------------------------------')
        print('dictionary: '+ uvar + '  keys:' + str(len(sp)))
        globals()[uvar]
        for ind in sp.index:
            try: 
                path = p['path'][ind].replace(', ','][')
                if p['type'][ind] == 'str': 
                    value = "'"+ str(p[site][ind]) + "'"
                else: value = str(p[site][ind])
                #exec(p['var'][ind] + p['path'][ind].replace(',','][').replace(' ', '') + ' = ' + str(p[site][ind]))
                line = uvar + path + ' = ' + value
                print(line)
                exec(line)
            #except: print('ERROR: '+ uvar + path + ' = ' + value + ' NOT UPDATED' )
            except Exception: traceback.print_exc()




""" 

### ------- values

wpara ={
    'swe':{
        'infolder': 'susiPath',
        'infile_d':'03-Korpis_weather.csv',
        'start_yr': 2000, 'end_yr': 2022, 
        'description': 'Somewhere in Sweden',
        'lat': 65.00, 'lon': 25.00
        }
}

cpara = {'dt': 86400.0,
        'flow' : { # flow field
                    'zmeas': 2.0,
                    'zground': 0.5,
                    'zo_ground': 0.01
                    },
        'interc': { # interception
                    'wmax': 0.5, # storage capacity for rain (mm/LAI)
                    'wmaxsnow': 4.0, # storage capacity for snow (mm/LAI),
                    },
        'snow': {
                # degree-day snow model
                'kmelt': 2.8934e-05, # melt coefficient in open (mm/s)
                'kfreeze': 5.79e-6, # freezing coefficient (mm/s)
                'r': 0.05 # maximum fraction of liquid in snow (-)
                },

        'physpara': {
                    # canopy conductance
                    'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                    'g1_conif': 2.1, # stomatal parameter, conifers
                    'g1_decid': 3.5, # stomatal parameter, deciduous
                    'q50': 50.0, # light response parameter (Wm-2)
                    'kp': 0.6, # light attenuation parameter (-)
                    'rw': 0.20, # critical value for REW (-),
                    'rwmin': 0.02, # minimum relative conductance (-)
                    # soil evaporation
                    'gsoil': 1e-2 # soil surface conductance if soil is fully wet (m/s)
                    },
        'phenopara': {
                    #seasonal cycle of physiology: smax [degC], tau[d], xo[degC],fmin[-](residual photocapasity)
                    'smax': 18.5, # degC
                    'tau': 13.0, # days
                    'xo': -4.0, # degC
                    'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                    },
        'state': {
                    'lai_conif': 3.0, # conifer 1-sided LAI (m2 m-2)
                    'lai_decid_max': 0.01, # maximum annual deciduous 1-sided LAI (m2 m-2): 
                    'hc': 16.0, # canopy height (m)
                    'cf': 0.7, # canopy closure fraction (-)
                    #initial state of canopy storage [mm] and snow water equivalent [mm]
                    'w': 0.0, # canopy storage mm
                    'swe': 0.0, # snow water equivalent mm
                    }
        }
org_para = {
        'org_depth': 0.04, # depth of organic top layer (m)
        'org_poros': 0.9, # porosity (-)
        'org_fc': 0.3, # field capacity (-)
        'org_rw': 0.24, # critical vol. moisture content (-) for decreasing phase in Ef
        'pond_storage_max': 0.01, # max ponding allowed (m)
        #initial states
        'org_sat': 1.0, # organic top layer saturation ratio (-)
        'pond_storage': 0.0 # pond storage
        }
    
# Hannun parametrit
#------------ Soil and stand parameters ----------------------------------
spara ={
    'swe':{
        'sitename': 'susirun',
        'species': 'Spruce', 'sfc':'sfc', 'sfc_specification': 1,
        'hdom':'hdomSim', 'vol':'volSim', 'age':'ageSim', 'smc': 'Peatland',
        'nLyrs':30, 'dzLyr': 0.05, 'L': 'sarkaSim', 'n':'n', 
        'ditch depth west': [-0.3, -0.9],   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
        'ditch depth east': [-0.3, -0.9],
        'ditch depth 20y west': [-0.3, -0.9],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'ditch depth 20y east': [-0.3, -0.9],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'scenario name': ['Control', 'DNM90'], #kasvunlisaykset
        'initial h': -0.2, 'slope': 3.0, 
        'peat type':['A','A','A','A','A','A','A','A'], 
        'peat type bottom':['A'],'anisotropy':10.,
        'vonP': True,
        'vonP top':  [2,2,2,2,4,5,5,5], 
        'vonP bottom': 5,
        'bd top':None, 'bd bottom': 0.16,
        'peatN':1.2, 'peatP':0.12, 'peatK':0.07,                               #peat nutrient contents in gravimetric %
        'enable_peattop': True, 'enable_peatmiddle': False,
        'enable_peatbottom': False,
        'depoN': 4.0, 'depoP':0.1, 'depoK':1.0,
        'fertilization': {
                'application year': 2201,
                'N':{'dose': 0.0, 'decay_k': 0.5, 'eff': 1.0},                              # fertilization dose in kg ha-1, decay_k in yr-1
                'P':{'dose': 45.0, 'decay_k': 0.2, 'eff': 1.0},
                'K':{'dose': 100.0, 'decay_k': 0.3, 'eff': 1.0},
                'pH_increment':1.0},  
    }
}

site = '03_Korpis'
param_file = "susi_2022\inputs\sweden\parameters.xlsx"


getVariablesInFile(param_file, site)
getDicInFile(param_file, site)
#spara['swe']['ditch depth 20y east']

#wpara['swe']['infile_d']

 """


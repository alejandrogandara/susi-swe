#%% imports
import pandas as pd
from datetime import datetime
from susi_2022.susi.figures import *
from netCDF4 import Dataset
from glob import glob
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress


def mk(path):
    if not os.path.exists(path): os.mkdir(path)
    return path

#%% initial parameters
scen = 0
#ind = 12

fs = 15  #graph font
facecolor = '#f2f5eb'
midLineColor = 'lightgray'
midLineStyle = 'dashed'
midLineWidth = 1
r2lim = 0.5
plim = 0.05

os.chdir('o:\\projects\\forestProductivity\\git\\susi-swe')

# Site description, extra fields
site_description_file = 'inputs/sweden/site_type.csv'
paramFile = pd.read_excel("inputs\sweden\parameters.xlsx")
sites = paramFile.columns[paramFile.columns.get_loc('source')+ 1:]

# Get the general parameters
working_folder = 'outputs'  #where the nc data is located
output_folder_graph = mk(f'{working_folder}/graphs')
files = glob(working_folder+'/*.nc')  # it will read all the nc files in folder
tfiles = pd.DataFrame({'files': files})

# Parameter file (use the exact one used to run the simulation)
tsites = pd.DataFrame({'sites':sites})
site_list = tfiles.join(tsites)



#%%  Individual Graphs
# if True == True:
#     ind= 0
#     wsite = sites[ind]

def stand_graph (level_1, level_2, fileName, output_folder_graph=output_folder_graph):
    for ind, wsite in enumerate(site_list.sites):
        #output folder
        output_folder = mk(f'{output_folder_graph}/{wsite}/')
        print(f'\nsite: {wsite}')

        
        try: 
            ncfFile = files[ind]

            ncf = Dataset(ncfFile, mode='r')  # Read the netCDF file
            print(ncfFile)
        except:
            print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
            continue

        basalarea = ncf.groups[level_1].variables[level_2][:]
        units = ncf.groups[level_1].variables[level_2].units
        basalarea_avg = np.mean(basalarea, axis=2)

        # Plot the basal area for each scenario
        for i in range(basalarea_avg.shape[0]):
            plt.plot(range(basalarea_avg.shape[1]), basalarea_avg[i,:], label=f"Scenario {i+1}")

        ns = ncf.dimensions['nscens'].size
        ny = ncf.dimensions['nyrs'].size
        sinfo = f' site: {wsite}, years: {ny}, scen = {ns}'
        #plt.annotate(sinfo, xy=(x, y), xycoords='axes fraction', fontsize=9, ha='left', va='top')
         
        plt.xlabel("Year")
        plt.ylabel( f"{units}")
        plt.legend()
        plt.title(f'{wsite} {level_2} -- {ny} years , {ns} scen')
        
        if fileName != False: plt.savefig(output_folder + '/' + fileName, bbox_inches='tight')
        plt.show()

        ncf.close()
#%%
#stand_graph('stand','volumegrowth', 'volumegrowth.png')
stand_graph('stand','volume', 'volume.png')
#stand_graph('stand','volumegrowth', 'volumegrowth.png')
#stand_graph('strip','dwt', False)



#%%









#####   TESTING ZONE









#%%
from susi_2022.susi.figures import *
from netCDF4 import Dataset 

for ind, wsite in enumerate(site_list.sites):
    #output folder
    output_folder = mk(f'{output_folder_graph}/{wsite}/')
    print(f'\nsite: {wsite}')
    
    try: 
        ncfFile = files[ind]
        compare_scens(ncfFile)
        
        plt.savefig(output_folder + '/' + wsite + '_compare_scens.png', bbox_inches='tight')
        plt.show()

    except:
        print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
        continue


#scen = 1
# hydrology(ff, scen)
# stand(ff, scen)
# mass(ff, scen)
# carbon(ff, scen)
# compare_1(ff, [0,1])

    

#%%  TEST INIDIVIDUALS 
#run = False

if True == True:
    ind= 0
    wsite = sites[ind]

#for ind, wsite in enumerate(site_list.sites):
    #output folder
    output_folder = mk(f'{working_folder}/graphs/{wsite}/')
    print(f'\nsite: {wsite}')


    # Read the netCDF file
    try: 
        ncfFile = files[ind]

        nc = Dataset(ncfFile, mode='r')
        print(ncfFile)
    except:
        print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
#        continue
#    basalarea = ncf.groups['stand'].variables['basalarea'][:]
    nc
    ns = nc.dimensions['nscens'].size
    print(f'scenarios {ns}')
    
#%%
ncf 
np.mean(ncf['strip']['dwtyr_latesummer'][1,:, :], axis = 1)


# %%

basalarea = ncf['stand']['basalarea'][:]

time = range(basalarea.shape[:])

# Plot the basal area for each scenario
for i in range(basalarea.shape[0]):
    plt.plot(time, basalarea[i,:,:].flatten(), label=f"Scenario {i+1}")

# Add labels and legend
plt.xlabel('Time (years)')
plt.ylabel('Basal Area (m$^2$/ha)')
plt.title('Basal Area over Time for Multiple Scenarios')
plt.legend()

# Show the plot
plt.show()


# %%


# %%
ncfFile = files[0]

nc_file = Dataset(ncfFile, mode='r')
print(ncfFile)
        
#%%
#nscens_var = nc_file.groups['scen'].variables['nscens']
#nscens_names = nscens_var.long_name.split(', ')
# %%
ncf.close()


#%%

df = pd.read_excel("inputs\sweden\parameters.xlsx", sheet_name='parameters', header=None)
df[1]
#df = df[(df[2] == 'pah') & (df[3] == "['netcdf']")]


#df = df.iloc[:, 9:]


#df = df.T.reset_index()
#df.columns = ['ncFile', 'site']


print(df)

# %%

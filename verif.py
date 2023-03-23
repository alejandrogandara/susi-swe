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
#ind = 5

fs = 15  #graph font
facecolor = '#f2f5eb'
os.chdir('o:\\projects\\forestProductivity\\git\\susi-swe')


gwl_observations = f'O:/projects/forestProductivity/01_data_acquisition/Ulf/SecondData/GVN_Man_ALL tidy by Alejandro Gandara.csv'
# Filter by the "Transekt" = 'mean', exclude transekts 1,2 and 3
obs_GWL = pd.read_csv(gwl_observations, encoding='latin1', sep=';')
obs_GWL = obs_GWL.loc[(obs_GWL.Transekt == 'mean')][['site','date', 'wt']]



# general parameters

working_folder = 'outputs'  #where the nc data is located
output_folder_graph = mk(f'{working_folder}/graphs')
#outpara = working_folder
files = glob(working_folder+'/*.nc')  # it will read all the nc files in folder
tfiles = pd.DataFrame({'files': files})


pFile = pd.read_excel("inputs\sweden\parameters.xlsx")  #calls the parameter file used to run susi, to collect initial date and site names
sites = pFile.columns[9:]

tsites = pd.DataFrame({'sites':sites})
tfiles.join(tsites)
print(tfiles.join(tsites))

# Calculated variables


#wsite = sites[ind]

#%%  Individual Graphs
#run = False
#if run = True:
for ind, wsite in enumerate(sites):
    #output folder
    output_folder = mk(f'{working_folder}/graphs/{wsite}/')

    print(f'site: {wsite}')

    #pFile[wsite]
    start_date = pFile.loc[pFile['key']=='start_date', [wsite] ].values[0]
    end_date = pFile.loc[pFile['key']=='end_date', [wsite] ].values[0]
    days = end_date - start_date
    years = end_date[0].year - start_date[0].year
    print(f"{years} years, and {days[0].days} days, from: {start_date[0].strftime('%Y-%m-%d')} to: {end_date[0].strftime('%Y-%m-%d')}")

    start_date = start_date[0].strftime('%Y-%m-%d')
    end_date = end_date[0].strftime('%Y-%m-%d')
    #ncf = Dataset('outputs\\05_StraRed_susi.nc', mode='r')
    try: 
        ncfFile = files[ind]

        ncf = Dataset(ncfFile, mode='r')
        print(ncfFile)
    except:
        print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
        continue

    # Filtered by the "Transekt" = 'mean',  transekts 1,2 and 3 area excluded
    obs_GWL = obs_GWL.loc[(obs_GWL.site == wsite)][['site','date', 'wt']]
    obs_GWL['date'] =  pd.to_datetime(obs_GWL['date'])
    obs_GWL['relative_day'] = obs_GWL['date'].apply(lambda x: (x - np.datetime64(start_date)).days)
    obs_GWL.set_index('relative_day', drop=False, inplace=True)
    
    print(len(obs_GWL))

    #Get water table
    #daily water table
    wt = np.mean(ncf['strip']['dwt'][scen,:, :], axis = 1)  

    #observed water table
    df_wt = pd.DataFrame({'wt':wt})  
    df_wt['relative_day'] = df_wt.index

    #Observed and modeled water table
    wt_comp = df_wt.join(obs_GWL, lsuffix='_est')  
    wt_comp['wt'] = pd.to_numeric(wt_comp['wt'], errors='coerce')
    wt_comp['date_est'] = pd.to_datetime(start_date) + pd.to_timedelta(wt_comp['relative_day_est'], unit='d')
    wt_comp = wt_comp.set_index('date_est')


    #    df_wt = pd.DataFrame({'wt':wt})
    #    df_wt['relative_day'] = df_wt.index
    #    wt_comp = df_wt.join(obs_GWL, lsuffix='_est')


    #    wt_comp['wt'] = pd.to_numeric(wt_comp['wt'], errors='coerce')
    #    wt_comp['date_est'] = pd.to_datetime(start_date) + pd.to_timedelta(wt_comp['relative_day_est'], unit='d')
    #    wt_comp = wt_comp.set_index('date_est')


    # Daily water table
    wt_fig = plt.figure(num='Water table', figsize=(15,3)) 
    if len(obs_GWL) > 0: 
        plt.scatter(wt_comp.index, wt_comp['wt'], label='wt_observed', color='darkblue')
    plt.plot(wt_comp['wt_est'], label='wt_estimated', color='green')
    #plt.xlabel(wt_comp['date_est'])
    plt.legend(loc=2)
    plt.title('Daily water table ' + wsite, fontsize = 15)
    plt.ylabel('WT m')
    plt.savefig(output_folder +  wsite + '_wt.png', bbox_inches='tight')
    #plt.show()
    plt.close()
    
    # Daily water table zoomed
    if len(obs_GWL) > 0:
        wt_fig_section = plt.figure(num='Water table zoom', figsize=(10,3)) 
        if len(obs_GWL) > 0: 
            tail = 200
            wt_comp_tail= wt_comp.loc[(wt_comp['relative_day_est'] >= (obs_GWL.relative_day.min()-tail)) & (wt_comp['relative_day_est'] <= (obs_GWL.relative_day.max()+tail))]
            plt.scatter(wt_comp_tail.index, wt_comp_tail['wt'], label='wt_observed', color='darkblue')
        plt.plot(wt_comp_tail['wt_est'], label='wt_estimated', color='green')
        #plt.xlabel(wt_comp['date_est'])
        plt.legend(loc=2)
        plt.title('Daily water table ' + wsite + ' (zoomed section)', fontsize = 15)
        plt.ylabel('WT m')
        plt.savefig(output_folder +  wsite + '_wt_section.png', bbox_inches='tight')
        #plt.show()
        plt.close()


        # read data from a CSV file
        data = wt_comp.dropna(subset=['wt'])
        # calculate linear regression line
        try: 
            slope, intercept, r_value, p_value, std_err = linregress(data['wt'], data['wt_est'])
            line = slope * data['wt'] + intercept


            xyfig = plt.figure(num='xy comp') 
            # plot observed vs. modeled water table with linear regression line
            plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
            plt.plot(data['wt'], line, color='red')
            plt.xlabel('Observed Water Table')
            plt.ylabel('Modeled Water Table')
            plt.title(wsite + '\nObserved vs. Modeled Water Table')

            # calculate and print R-squared value as a measure of model fit
            r_squared = r_value ** 2
            #stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
            stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'

            plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

            # show plot and save
            #plt.show()
            xyfig.savefig(output_folder + wsite + '_wt_observed vs modeled.png', bbox_inches='tight')

            data.to_csv(output_folder + wsite + '_wt_observed_modeled.csv', sep=';')

            #  
            lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
                            'stderr':std_err}, index=[0])
            lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')
            plt.close()
        except: print("no linear model")
    ncf.close()

#%% all togethere
# Initialize an empty list to store the dataframes
df_list = []

# Loop through all subdirectories and files with the "_wt_lm.csv" extension
for root, dirs, files in os.walk(working_folder):
    for file in files:
        if file.endswith('wt_observed_modeled.csv'):
            # Read the csv file into a pandas dataframe
            df = pd.read_csv(os.path.join(root, file), sep=',|;', engine='python')
            # Append the dataframe to the list
            df_list.append(df)

# Concatenate all dataframes in the list into a single dataframe
consolidated_df = pd.concat(df_list)
consolidated_df = consolidated_df.dropna()

# Read the description.csv file into a pandas dataframe
site_description = pd.read_csv('inputs/sweden/site_type.csv', sep=',|;', encoding='latin1', engine='python')

# Merge the consolidated dataframe with the description dataframe on the "site" column
data = pd.merge(consolidated_df, site_description, on='site')



#import matplotlib.pyplot as plt

# calculate linear regression line
slope, intercept, r_value, p_value, std_err = linregress(data['wt'], data['wt_est'])
line = slope * data['wt'] + intercept


#%%  xy graph all sites
xyfig = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
plt.plot(data['wt'], line, color='red')
plt.xlabel('Observed Water Table')
plt.ylabel('Modeled Water Table')
plt.title('All sites\nObserved vs. Modeled Water Table')

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

# show plot and save
plt.show()
xyfig.savefig(output_folder_graph + '/all_water_table_observed_modeled.png', bbox_inches='tight')

data.to_csv(output_folder_graph + '/all_water_table_observed_modeled.csv', sep=';')

# #  
# lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
#                    'stderr':std_err}, index=[0])
# lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')






#%%  xy graph all sitesn (soil and species)
xyfig2 = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
plt.plot(data['wt'], line, color='red')
plt.xlabel('Observed Water Table')
plt.ylabel('Modeled Water Table')
plt.title('All sites\nObserved vs. Modeled Water Table')

# Group the merged dataframe by "soil_type" and "specie"
grouped_df = data.groupby(['soil_type', 'species'])

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o')
plt.legend(loc='lower right')
    

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

# show plot and save
plt.show()
xyfig2.savefig(output_folder_graph + '/all_wt_observed vs modeled_soil_spec.png', bbox_inches='tight')

#data.to_csv(output_folder +  wsite + '_wt_observed_modeled.csv', sep=';')

#  
#lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
#                   'stderr':std_err}, index=[0])
#lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')

#%%  xy graph all sites (soil type)
xyfig3 = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
plt.plot(data['wt'], line, color='red')
plt.xlabel('Observed Water Table')
plt.ylabel('Modeled Water Table')
plt.title('All sites\nObserved vs. Modeled Water Table')

grouped_df = data.groupby(['soil_type'])

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o')
plt.legend(loc='lower right')

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

# show plot and save
plt.show()
xyfig3.savefig(output_folder_graph + '/all_wt_observed vs modeled_soil.png', bbox_inches='tight')

#%%  xy graph all sites (species)
xyfig4 = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
plt.plot(data['wt'], line, color='red')
plt.xlabel('Observed Water Table')
plt.ylabel('Modeled Water Table')
plt.title('All sites\nObserved vs. Modeled Water Table')

grouped_df = data.groupby(['species'])

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o')
plt.legend(loc='lower right')

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

# show plot and save
plt.show()
xyfig4.savefig(output_folder_graph + '/all_wt_observed vs modeled_species.png', bbox_inches='tight')


#%% Gropued graph p-value and R2
# Initialize an empty list to store the dataframes
df_list = []

# Loop through all subdirectories and files with the "_wt_lm.csv" extension
for root, dirs, files in os.walk(working_folder):
    for file in files:
        if file.endswith('_wt_lm.csv'):
            # Read the csv file into a pandas dataframe
            df = pd.read_csv(os.path.join(root, file), sep=';')
            # Append the dataframe to the list
            df_list.append(df)

# Concatenate all dataframes in the list into a single dataframe
consolidated_df = pd.concat(df_list)
site_description = pd.read_csv('inputs/sweden/site_type.csv', sep=';', encoding='latin1')

# Merge the consolidated dataframe with the description dataframe on the "site" column
merged_df = pd.merge(consolidated_df, site_description, on='site')

# Group the merged dataframe by "soil_type" and "specie"
grouped_df = merged_df.groupby(['soil_type', 'species'])
fig, ax = plt.subplots()

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o')

# Add the legend and axis labels
ax.legend()
ax.set_xlabel('pvalue')
ax.set_ylabel('r_squared')

# Display the plot
plt.savefig(output_folder_graph + '/all_gstats.png', bbox_inches='tight')
plt.show()

# Grouped by SoilType
grouped_df = merged_df.groupby(['soil_type'])

# Initialize the plot
fig, ax = plt.subplots()

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o')

# Add the legend and axis labels
ax.legend()
ax.set_xlabel('pvalue')
ax.set_ylabel('r_squared')

# Display the plot
plt.savefig(output_folder_graph + '/all_gstats2.png', bbox_inches='tight')
plt.show()

# graph
grouped_df = merged_df.groupby(['species'])
fig, ax = plt.subplots()

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o')
ax.legend()
ax.set_xlabel('pvalue')
ax.set_ylabel('r_squared')


# Display the plot
plt.savefig(output_folder_graph + '/all_gstats3.png', bbox_inches='tight')
plt.show()

# %%

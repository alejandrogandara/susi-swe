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

# Water table observation file
gwl_observations = f'O:/projects/forestProductivity/01_data_acquisition/Ulf/SecondData/GVN_Man_ALL tidy by Alejandro Gandara.csv'
# Filter by the "Transekt" = 'mean', exclude transekts 1,2 and 3
obs_GWL_all = pd.read_csv(gwl_observations, encoding='latin1', sep=';')
obs_GWL_all = obs_GWL_all.loc[(obs_GWL_all.Transekt == 'mean')][['site','date', 'wt']]



# Get the general parameters

working_folder = 'outputs'  #where the nc data is located
output_folder_graph = mk(f'{working_folder}/graphs')
files = glob(working_folder+'/*.nc')  # it will read all the nc files in folder
tfiles = pd.DataFrame({'files': files})

# Parameter file (use the exact one used to run the simulation)
pFile = pd.read_excel("inputs\sweden\parameters.xlsx")
sites = pFile.columns[9:]

tsites = pd.DataFrame({'sites':sites})
site_list = tfiles.join(tsites)
print(tfiles.join(tsites))

# Calculated variables
print(tfiles.join(tsites))

#wsite = sites[ind]


#%% temporal no funciona
# paramFile = pd.read_excel("inputs\sweden\parameters.xlsx")

# source_col_index = paramFile.columns.get_loc('source')
# column_names = list(paramFile.columns[source_col_index+1:])
# column_names

# dominant_age = float(paramFile.loc[paramFile['var'] == "age_dominant"][site].values[0])


#%%  Individual Graphs
#run = False

# if True == True:
#     ind= 1
#     wsite = sites[ind]

for ind, wsite in enumerate(site_list.sites):
    #output folder
    output_folder = mk(f'{working_folder}/graphs/{wsite}/')
    print(f'\nsite: {wsite}')

    #pFile[wsite]
    # Adjust relative dates in the nc file using the initial date.
    start_date = pFile.loc[pFile['key']=='start_date', [wsite] ].values[0]
    end_date = pFile.loc[pFile['key']=='end_date', [wsite] ].values[0]
    days = end_date - start_date
    years = end_date[0].year - start_date[0].year
    print(f"{years} years, and {days[0].days} days, from: {start_date[0].strftime('%Y-%m-%d')} to: {end_date[0].strftime('%Y-%m-%d')}")

    start_date = start_date[0].strftime('%Y-%m-%d')
    end_date = end_date[0].strftime('%Y-%m-%d')

    # Read the netCDF file
    try: 
        ncfFile = files[ind]

        ncf = Dataset(ncfFile, mode='r')
        print(ncfFile)
    except:
        print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
        continue

    # Filtered the Observed Water Table (m) by the "Transekt" = 'mean',  transekts 1,2 and 3 are excluded
    obs_GWL = obs_GWL_all.loc[(obs_GWL_all.site == wsite)][['site','date', 'wt']]
    obs_GWL['date'] =  pd.to_datetime(obs_GWL['date'])
    obs_GWL['relative_day'] = obs_GWL['date'].apply(lambda x: (x - np.datetime64(start_date)).days)
    obs_GWL.set_index('relative_day', drop=False, inplace=True)
    
    print(f'water table valid observations: {len(obs_GWL)}')

    # Get water table from the netCDF file
    # Daily water table
    wt = np.mean(ncf['strip']['dwt'][scen,:, :], axis = 1)  

    # Modeled water table
    df_wt = pd.DataFrame({'wt':wt})  
    df_wt['relative_day'] = df_wt.index

    #Observed and modeled water table
    wt_comp = df_wt.join(obs_GWL, lsuffix='_est')  
    wt_comp['wt'] = pd.to_numeric(wt_comp['wt'], errors='coerce')
    wt_comp['date_est'] = pd.to_datetime(start_date) + pd.to_timedelta(wt_comp['relative_day_est'], unit='d')
    wt_comp = wt_comp.set_index('date_est')

    # Daily water table Graph  --------------------------------------------
    wt_fig = plt.figure(num='Water table', figsize=(15,3)) 
    if len(obs_GWL) > 0: 
        plt.scatter(wt_comp.index, wt_comp['wt'], label='wt_observed', color='darkblue')
    plt.plot(wt_comp['wt_est'], label='wt_estimated', color='green')
    #plt.xlabel(wt_comp['date_est'])
    #plt.axhline(min(wt_comp['wt']))
    
    # Reference mean lines
    wtObsMean = np.mean(wt_comp['wt'])
    wtEstMean = np.mean(wt_comp['wt_est'])
    
    plt.axhline(wtObsMean, color='lightblue', linestyle=midLineStyle, linewidth=midLineWidth, label= '_wt_observed_mean')
    plt.axhline(wtEstMean, color='palegreen', linestyle=midLineStyle, linewidth=midLineWidth, label= '_wt_estimated_mean')
    
    plt.legend(loc=2)
    plt.title('Daily water table ' + wsite, fontsize = 15)
    plt.ylabel('WT m')
    plt.savefig(output_folder +  wsite + '_wt.png', bbox_inches='tight')
    plt.show()
    #plt.close()
    
    # Daily water table zoomed graph -----------------------------------
    if len(obs_GWL) > 0:
        wt_fig_section = plt.figure(num='Water table zoom', figsize=(10,3)) 
        if len(obs_GWL) > 0: 
            tail = 200
            wt_comp_tail= wt_comp.loc[(wt_comp['relative_day_est'] >= (obs_GWL.relative_day.min()-tail)) & (wt_comp['relative_day_est'] <= (obs_GWL.relative_day.max()+tail))]
            plt.scatter(wt_comp_tail.index, wt_comp_tail['wt'], label='wt_observed', color='darkblue')
        plt.plot(wt_comp_tail['wt_est'], label='wt_estimated', color='green')

        plt.axhline(wtObsMean, color='lightblue', linestyle=midLineStyle, linewidth=midLineWidth, label= '_wt_observed_mean')
        plt.axhline(wtEstMean, color='palegreen', linestyle=midLineStyle, linewidth=midLineWidth, label= '_wt_estimated_mean')
        
        #plt.xlabel(wt_comp['date_est'])
        plt.legend(loc=2)
        plt.title('Daily water table ' + wsite + ' (zoomed section)', fontsize = 15)
        plt.ylabel('WT m')
        plt.savefig(output_folder +  wsite + '_wt_section.png', bbox_inches='tight')
        plt.show()
        #plt.close()

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
            plt.xlabel('Observed Water Table (m)')
            plt.ylabel('Modeled Water Table (m)')
            plt.title(wsite + '\nObserved vs. Modeled Water Table')

            # calculate and print R-squared value as a measure of model fit
            r_squared = r_value ** 2
            #stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
            stats_text = f'R\u00b2: {r_squared:.2f}  p-value: {p_value:.2f}'
            plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

            #Midline
            xymin= min([min(data['wt']), min(data['wt_est'])])
            plt.plot([xymin, 0], [xymin, 0], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
            plt.xlim([xymin, 0])
            plt.ylim([xymin, 0])

            # show plot and save
            #plt.show()
            xyfig.savefig(output_folder + wsite + '_wt_observed vs modeled.png', bbox_inches='tight')

            data.to_csv(output_folder + wsite + '_wt_observed_modeled.csv', sep=';')

            #  
            lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
                            'stderr':std_err}, index=[0])
            lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')
            plt.show()
            #plt.close()
        except: print("no linear model")
    ncf.close()

#%% all together
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
# calculate residuals
data['residuals'] = data['wt_est'] - data['wt']
rmse = np.sqrt(np.mean(data['residuals']**2))

#import matplotlib.pyplot as plt
#%%
# calculate linear regression line
slope, intercept, r_value, p_value, std_err = linregress(data['wt'], data['wt_est'])
line = slope * data['wt'] + intercept

xyfig = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
plt.plot(data['wt'], line, color='red')
plt.xlabel('Observed Water Table (m)')
plt.ylabel('Modeled Water Table (m)')
plt.title('All sites\nObserved vs. Modeled Water Table')

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
#stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')


# show plot and save
plt.show()
xyfig.savefig(output_folder_graph + '/all_water_table_observed_modeled.png', bbox_inches='tight')

data.to_csv(output_folder_graph + '/wtlm_AllSites.csv', sep=';')

# #  
# lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
#                    'stderr':std_err}, index=[0])
# lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')




#%%  xy graph all sites (species)
wtlm_species = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
grouped_df = data.groupby(['species'])
x = 0.05
y = 1
names = ''
# Loop through each group and fit a linear regression model
for name, group in grouped_df:
    # Fit the model
    slope, intercept, r_value, p_value, std_err = linregress(group['wt'], group['wt_est'])
    line = slope * group['wt'] + intercept
    rmse = np.sqrt(np.mean(group['residuals']**2))

    # Plot the data and the fitted line
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o', alpha=0.5)
    plt.plot(group['wt'], line, label="_name")

    # Annotate the plot with the regression statistics
    y = y - 0.04
    stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}  {name}'
    plt.annotate(stats_text, xy=(x, y), xycoords='axes fraction', fontsize=9, ha='left', va='top')
    names = name + ', ' + names
    # Set the axis labels and title

#Midline
xymin= min([min(group['wt']), min(group['wt_est'])])
plt.plot([xymin, 0], [xymin, 0], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
plt.xlim([xymin, 0])
plt.ylim([xymin, 0])

#Titles
plt.xlabel('Observed Water Table (m)')
plt.ylabel('Modeled Water Table (m)')
plt.title(f'{names}\nObserved vs. Modeled Water Table')
plt.legend(loc='lower right')
plt.show()

wtlm_species.savefig(output_folder_graph + '/wtlm_species.png', bbox_inches='tight')


#%%  xy graph all sites (SOIL TYPE)
wtlm_soilType = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
grouped_df = data.groupby(['soil_type'])
x = 0.05
y = 1
names = ''
# Loop through each group and fit a linear regression model
for name, group in grouped_df:
    # Fit the model
    slope, intercept, r_value, p_value, std_err = linregress(group['wt'], group['wt_est'])
    line = slope * group['wt'] + intercept
    rmse = np.sqrt(np.mean(group['residuals']**2))

    # Plot the data and the fitted line
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o', alpha=0.5)
    plt.plot(group['wt'], line, label="_name")

    # Annotate the plot with the regression statistics
    y = y - 0.04
    stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}  {name}'
    
    plt.annotate(stats_text, xy=(x, y), xycoords='axes fraction', fontsize=9, ha='left', va='top')
    names = name + ', ' + names
    # Set the axis labels and title
    
#Midline
xymin= min([min(group['wt']), min(group['wt_est'])])
plt.plot([xymin, 0], [xymin, 0], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
plt.xlim([xymin, 0])
plt.ylim([xymin, 0])

plt.xlabel('Observed Water Table (m)')
plt.ylabel('Modeled Water Table (m)')
plt.title(f'{names}\nObserved vs. Modeled Water Table')
plt.legend(loc='lower right')
plt.show()

wtlm_soilType.savefig(output_folder_graph + '/wtlm_soylType.png', bbox_inches='tight')


#%%  xy graph all sites (SOIL TYPE)
wtlm_location = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
#plt.scatter(data['wt'], data['wt_est'], alpha=0.5)
grouped_df = data.groupby(['side'])
x = 0.05
y = 1
names = ''
# Loop through each group and fit a linear regression model
for name, group in grouped_df:
    # Fit the model
    slope, intercept, r_value, p_value, std_err = linregress(group['wt'], group['wt_est'])
    line = slope * group['wt'] + intercept
    rmse = np.sqrt(np.mean(group['residuals']**2))

    # Plot the data and the fitted line
    plt.scatter(group['wt'], group['wt_est'], label=name, marker='o', alpha=0.5)
    plt.plot(group['wt'], line, label="_name")

    # Annotate the plot with the regression statistics
    y = y - 0.04
    stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}  {name}'
    plt.annotate(stats_text, xy=(x, y), xycoords='axes fraction', fontsize=9, ha='left', va='top')
    names = name + ', ' + names
    # Set the axis labels and title

#Midline
xymin= min([min(group['wt']), min(group['wt_est'])])
plt.plot([xymin, 0], [xymin, 0], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
plt.xlim([xymin, 0])
plt.ylim([xymin, 0])

plt.xlabel('Observed Water Table (m)')
plt.ylabel('Modeled Water Table (m)')
plt.title(f'{names}\nObserved vs. Modeled Water Table')
plt.legend(loc='lower right')
plt.show()

wtlm_location.savefig(output_folder_graph + '/wtlm_location.png', bbox_inches='tight')


#%% Gropued graph p-value and R2
# Initialize an empty list to store the dataframes
df_list = []



# Loop through all subdirectories and files with the "_wt_lm.csv" extension
for root, dirs, files in os.walk(working_folder):
    for file in files:
        if file.endswith('_wt_lm.csv'):
            df = pd.read_csv(os.path.join(root, file), sep=';')
            df_list.append(df)

# Concatenate all dataframes in the list into a single dataframe
consolidated_df = pd.concat(df_list)
site_description = pd.read_csv(site_description_file, sep=';', encoding='latin1')

# Merge the consolidated dataframe with the description dataframe on the "site" column
merged_df = pd.merge(consolidated_df, site_description, on='site')

# Grouped by SoilType
grouped_df = merged_df.groupby(['soil_type'])

# Initialize the plot
fig, ax = plt.subplots()

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o', alpha=0.5)
    for i in range(len(group)):
        ax.annotate('   ' + group['site'].iloc[i], 
                (group['pvalue'].iloc[i], group['r_squared'].iloc[i]), fontsize=8)

# reference line
ax.axvline(plim, color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
ax.axhline(r2lim, color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)

# Add the legend and axis labels
ax.legend(loc='upper right')
ax.set_xlabel('pvalue')
ax.set_ylabel('r_squared')
ax.set_title('Observed vs. Modeled Water Table')
plt.xlim([0, 0.2])
plt.ylim([0, 1])


# Display the plot
plt.savefig(output_folder_graph + '/stats_soilType.png', bbox_inches='tight')
plt.show()

#%%  PINE AND SPRUCE GRAPH
# graph
grouped_df = merged_df.groupby(['species'])
fig, ax = plt.subplots()

# Loop through each group and plot the data with a different shape
for name, group in grouped_df:
    ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o', alpha=0.5)
    for i in range(len(group)):
        ax.annotate('   ' + group['site'].iloc[i], 
            (group['pvalue'].iloc[i], group['r_squared'].iloc[i]), fontsize=8)

# reference line
ax.axvline(plim, color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
ax.axhline(r2lim, color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)


ax.legend(loc='upper right')
ax.set_xlabel('pvalue')
ax.set_ylabel('r_squared')
ax.set_title('Observed vs. Modeled Water Table')

plt.xlim([0, 0.2])
plt.ylim([0, 1])


plt.savefig(output_folder_graph + '/stats_species.png', bbox_inches='tight')
plt.show()



#%%  this is only to asjust the 
# for name, group in grouped_df:
#     ax.scatter(group['pvalue'], group['r_squared'], label=name, marker='o')
#     for i in range(len(group)):
#         ax.annotate(group['site'].iloc[i] + ' (' + str(group['lat'].iloc[i]) + ')', 
#                     (group['pvalue'].iloc[i], group['r_squared'].iloc[i]))




# # %%
# increment = [0, 0.5, 1]
# initial = -0.6

# y  = [round(initial * (i + 1),2) for i in increment]
# y2 = [-0.01] * 5
# y
# lab  = ['scen_' + str((i + 1)) for i in increment]

# # %%
# y2
# # %%
# lab
# # %%

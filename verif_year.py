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

#netCDF file settings
level_1 = 'stand'  #in nc file
level_2 = 'basalarea' #in nc file

# labels not related to the file
var_label = 'BA'  
var_label_long = 'Basal Area'

#graph fix params
fs = 15  
facecolor = '#f2f5eb'
midLineColor = 'lightgray'
midLineStyle = 'dashed'
midLineWidth = 1
r2lim = 0.5
plim = 0.05


os.chdir('o:\\projects\\forestProductivity\\git\\susi-swe')

working_folder = 'outputs'  #where the nc data is located
files = glob(working_folder+'/*.nc')  # it will read all the nc files in folder
tfiles = pd.DataFrame({'files': files})

output_folder_graph = mk(f'{working_folder}/graphs')

# Site description, extra fields from excel file
site_description_file = 'inputs/sweden/site_type.csv'
paramFile = pd.read_excel("inputs/sweden/parameters.xlsx")

# Observed variable to compare with the nc file, in excel file
bai_observations = pd.read_excel(f'O:/projects/forestProductivity/01_data_acquisition/Ulf/BAI/BAI.xlsx')

# Filter the useful columns, match yor field names in excel and double check 'var_label'.
df_var_obs_all = bai_observations[['site', 'year', var_label]]



# Get the general parameters


# Parameter file (use the exact one used to run the simulation)
pFile = pd.read_excel("inputs\sweden\parameters.xlsx")

sites = pFile.columns[pFile.columns.get_loc('source')+ 1:]

tsites = pd.DataFrame({'sites':sites})
site_list = tfiles.join(tsites)
print(tfiles.join(tsites))

# Calculated variables
print(tfiles.join(tsites))

#wsite = sites[ind]


#%%  Individual Graphs
#run = False

# if True == True:
#     ind= 0
#     wsite = sites[ind]

for ind, wsite in enumerate(site_list.sites):
    #output folder
    output_folder = mk(f'{working_folder}/graphs/{wsite}/')
    print(f'\nsite: {wsite}')

    #pFile[wsite]
    # Adjust relative dates in the nc file using the initial date.
    start_date = pFile.loc[pFile['key']=='start_date', [wsite] ].values[0]
    start_year = start_date[0].year
       
    end_date = pFile.loc[pFile['key']=='end_date', [wsite] ].values[0]
    end_year = end_date[0].year
    
    print(f"{wsite} {end_year - start_year} years, from: {start_year} to: {end_year}")

    # Read the netCDF file
    try: 
        ncfFile = files[ind]

        ncf = Dataset(ncfFile, mode='r')
        units = ncf.groups[level_1].variables[level_2].units
        print(ncfFile)
    except:
        print(f'ncf file not found or damaged {wsite} - {ncfFile}\n')
        continue

    # Filtered the Observed Water Table (m) by the "Transekt" = 'mean',  transekts 1,2 and 3 are excluded
    df_var_obs = df_var_obs_all.loc[(df_var_obs_all.site == wsite)][['site','year', var_label]]
#--    df_var_obs['year'] =  pd.to_datetime(df_var_obs['year'])
    df_var_obs['relative_year'] = df_var_obs['year'].apply(lambda x: (x - start_year))
    df_var_obs.set_index('relative_year', drop=False, inplace=True)
    
    print(f'estimated BA using tree cores: {len(df_var_obs)}')
    
    # Get Variable from the netCDF file
    # Yearly or daily values
    var = np.mean(ncf[level_1][level_2][scen,:, :], axis = 1)  

    # Modeled variable
    df_var_model = pd.DataFrame({var_label:var})  
    df_var_model['relative_year'] = df_var_model.index
    df_var_model.set_index('relative_year', drop=False, inplace=True)

    #Observed and modeled
    var_est = f'{var_label}_est'
    var_comp = df_var_model.join(df_var_obs, lsuffix='_est', how='outer')    
    var_comp[var_label] = pd.to_numeric(var_comp[var_label], errors='coerce')
    var_comp['date_est'] = start_year + var_comp.index
    #var_comp['date_est'] = pd.to_datetime(start_date) + pd.to_timedelta(var_comp['relative_year_est'], unit='d')
    var_comp = var_comp.set_index('date_est')
    
    # Daily water table Graph  --------------------------------------------
    wt_fig = plt.figure(num=level_2, figsize=(15,3)) 
    if len(df_var_obs) > 0: 
        plt.scatter(var_comp.index, var_comp[var_label], label='observed', color='darkblue')
    plt.plot(var_comp[var_est], label='estimated', color='green')
    #plt.xlabel(var_comp['date_est'])
    #plt.axhline(min(var_comp[var_label]))
    
    # Reference mean lines
    wtObsMean = np.mean(var_comp[var_label])
    wtEstMean = np.mean(var_comp[var_est])
    
    plt.axhline(wtObsMean, color='lightblue', linestyle=midLineStyle, linewidth=midLineWidth, label= '_observed_mean')
    plt.axhline(wtEstMean, color='palegreen', linestyle=midLineStyle, linewidth=midLineWidth, label= '_estimated_mean')
    
    plt.legend(loc=2)
    plt.title(level_2 + ' ' + wsite, fontsize = 15)
    plt.ylabel('WT m')
    plt.savefig(output_folder +  wsite + '_' + level_2 + '.png', bbox_inches='tight')
    plt.show()
    #plt.close()
    
    # Daily water table zoomed graph -----------------------------------
    if len(df_var_obs) > 0:
        wt_fig_section = plt.figure(num='zoomed', figsize=(10,3)) 
        if len(df_var_obs) > 0: 
            tail = 200
            var_comp_tail= var_comp.loc[(var_comp['relative_year_est'] >= (df_var_obs.relative_year.min()-tail)) & (var_comp['relative_year_est'] <= (df_var_obs.relative_year.max()+tail))]
            plt.scatter(var_comp_tail.index, var_comp_tail[var_label], label='observed', color='darkblue')
        plt.plot(var_comp_tail[var_est], label='estimated', color='green')

        plt.axhline(wtObsMean, color='lightblue', linestyle=midLineStyle, linewidth=midLineWidth, label= '_observed_mean')
        plt.axhline(wtEstMean, color='palegreen', linestyle=midLineStyle, linewidth=midLineWidth, label= '_estimated_mean')
        
        #plt.xlabel(var_comp['date_est'])
        plt.legend(loc=2)
        plt.title(wsite +' '+ units + ' (zoomed section)', fontsize = 15)
        plt.ylabel('WT m')
        plt.savefig(output_folder +  wsite + '_'+ level_2 +'_section.png', bbox_inches='tight')
        plt.show()
        #plt.close()

       # read data from a CSV file
        data_ = var_comp.dropna(subset=[var_label, var_est])
        data = data_[1:]
        
        #calculate residuals
        data['residuals'] = data[var_est] - data[var_label]
        rmse = np.sqrt(np.mean(data['residuals']**2))
        
        
        # calculate linear regression line
        try: 
            slope, intercept, r_value, p_value, std_err = linregress(data[var_label], data[var_est])
            line = slope * data[var_label] + intercept

            xyfig = plt.figure(num='xy comp') 
            # plot observed vs. modeled water table with linear regression line
            plt.scatter(data[var_label], data[var_est], alpha=0.5)
            plt.plot(data[var_label], line, color='red')
            plt.xlabel(f'Observed {units}')
            plt.ylabel(f'Modeled {units}')
            plt.title(wsite +' '+ var_label_long + '\nObserved vs. Modeled')

            # calculate and print R-squared value as a measure of model fit
            r_squared = r_value ** 2
            #stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
            stats_text = f'R\u00b2: {r_squared:.2f}  p-val: {p_value:.2f}, RMSE:{rmse:.2f}'
            plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')

            #Midline
            xymin= min([min(data[var_label]), min(data[var_est])])
            xymax= max([max(data[var_label]), max(data[var_est])])
            xymax = xymax * 1.01
            
            plt.plot([xymin, xymax], [xymin, xymax], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
            plt.xlim([xymin, xymax])
            plt.ylim([xymin, xymax])

            # show plot and save
            #plt.show()
            xyfig.savefig(output_folder + wsite + '_' + level_2 + '_observed vs modeled.png', bbox_inches='tight')

            data.to_csv(output_folder + wsite +  '_' + level_2 +'_observed_modeled.csv', sep=';')

            #  
            lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
                            'stderr':std_err}, index=[0])
            lm_restult.to_csv(output_folder +  wsite +'_'+ level_2 + '_lm.csv', sep=';')
            plt.show()
            #plt.close()
        except: print("no linear model")
    #ncf.close()

#%% all together
# Initialize an empty list to store the dataframes
df_list = []

# Loop through all subdirectories and files with the "_wt_lm.csv" extension
for root, dirs, files in os.walk(working_folder):
    for file in files:
        if file.endswith(f'{var_label}_observed_modeled.csv'):
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
data['residuals'] = data[var_est] - data[var_label]
rmse = np.sqrt(np.mean(data['residuals']**2))


#%%
#import matplotlib.pyplot as plt

# calculate linear regression line
slope, intercept, r_value, p_value, std_err = linregress(data[var_label], data[var_est])
line = slope * data[var_label] + intercept

xyfig = plt.figure(num='xy comp') 
# plot observed vs. modeled water table with linear regression line
plt.scatter(data[var_label], data[var_est], alpha=0.5)
plt.plot(data[var_label], line, color='red')
plt.xlabel(f'Observed {units}')
plt.ylabel(f'Modeled {units}')
plt.title(f'All sites, {var_label_long}\nObserved vs. Modeled')

# calculate and print R-squared value as a measure of model fit
r_squared = r_value ** 2
#stats_text = f'Slope: {slope:.2f}\nIntercept: {intercept:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
#stats_text = f'std error: {std_err:.2f}\nR-squared: {r_squared:.2f}\np-value: {p_value:.2f}'
stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}'
plt.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, ha='left', va='top')



#Midline
xymin= min([min(data[var_label]), min(data[var_est])])
xymax= max([max(data[var_label]), max(data[var_est])])
xymax = xymax * 1.01

plt.plot([xymin, xymax], [xymin, xymax], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
plt.xlim([xymin, xymax])
plt.ylim([xymin, xymax])


# show plot and save
plt.show()
xyfig.savefig(output_folder_graph + '/' + level_2 +'_all_observed_modeled.png', bbox_inches='tight')

data.to_csv(output_folder_graph + '/' + level_2 + '_AllSites.csv', sep=';')

# #  
# lm_restult = pd.DataFrame({'site':wsite, 'rvalue': r_value, 'r_squared' : r_squared,'pvalue': p_value,
#                    'stderr':std_err}, index=[0])
# lm_restult.to_csv(output_folder +  wsite + '_wt_lm.csv', sep=';')


#%%  xy graph all sites (species)
def obs_vs_model(groupby):
    wtlm_species = plt.figure(num='xy comp') 
    # plot observed vs. modeled water table with linear regression line
    #plt.scatter(data[var_label], data[var_est], alpha=0.5)
    grouped_df = data.groupby([groupby])
    x = 0.05
    y = 1
    names = ''
    # Loop through each group and fit a linear regression model
    for name, group in grouped_df:
        # Fit the model
        slope, intercept, r_value, p_value, std_err = linregress(group[var_label], group[var_est])
        line = slope * group[var_label] + intercept
        rmse = np.sqrt(np.mean(group['residuals']**2))

        # Plot the data and the fitted line
        plt.scatter(group[var_label], group[var_est], label=name, marker='o', alpha=0.5)
        plt.plot(group[var_label], line, label="_name")

        # Annotate the plot with the regression statistics
        y = y - 0.04
        stats_text = f'R\u00b2:{r_value ** 2:.2f},  p-val:{p_value:.2f}, RMSE:{rmse:.2f}  {name}'
        plt.annotate(stats_text, xy=(x, y), xycoords='axes fraction', fontsize=9, ha='left', va='top')
        names = name + ', ' + names
        # Set the axis labels and title

    #Midline
        xymin= min([min(data[var_label]), min(data[var_est])])
        xymax= max([max(data[var_label]), max(data[var_est])])
        xymax = xymax * 1.01

        plt.plot([xymin, xymax], [xymin, xymax], color=midLineColor, linestyle=midLineStyle, linewidth=midLineWidth)
        plt.xlim([xymin, xymax])
        plt.ylim([xymin, xymax])

    #Titles
    plt.xlabel(f'Observed {units}')
    plt.ylabel(f'Modeled {units}')
    plt.title(f'{names}\nObserved vs. Modeled {var_label_long}')
    plt.legend(loc='lower right')
    plt.show()

    wtlm_species.savefig(output_folder_graph + '/' + level_2 + '_'+ groupby + '.png', bbox_inches='tight')

#%%
obs_vs_model('species')
obs_vs_model('soil_type')
obs_vs_model('side')




#%% Gropued graph p-value and R2
# Initialize an empty list to store the dataframes
df_list = []

# Loop through all subdirectories and files with the "_wt_lm.csv" extension
for root, dirs, files in os.walk(working_folder):
    for file in files:
        if file.endswith(f'{var_label}_lm.csv'):
            df = pd.read_csv(os.path.join(root, file), sep=';')
            df_list.append(df)

# Concatenate all dataframes in the list into a single dataframe
consolidated_df = pd.concat(df_list)
site_description = pd.read_csv(site_description_file, sep=';', encoding='latin1')

# Merge the consolidated dataframe with the description dataframe on the "site" column
merged_df = pd.merge(consolidated_df, site_description, on='site')


#%%
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
ax.set_title(f'Observed vs. Modeled {var_label_long}')

plt.xlim([0, 0.2])
plt.ylim([0, 1])

# Display the plot
plt.savefig(output_folder_graph + '/'+ level_2 +'_stats_soilType.png', bbox_inches='tight')
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
ax.set_title(f'Observed vs. Modeled {var_label_long}')
plt.xlim([0, 0.2])
plt.ylim([0, 1])


plt.savefig(output_folder_graph + '/'+ level_2 +'_stats_species.png', bbox_inches='tight')
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

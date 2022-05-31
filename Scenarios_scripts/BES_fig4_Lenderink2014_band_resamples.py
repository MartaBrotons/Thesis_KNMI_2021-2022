# Python script to compute the un-resampled EC-Earth statistic + CMIP6 statistics (figure 4 from Lenderink 2014) (Figure 11 and 12)

# To make it work the output from the esmvaltool scenarios must be uploaded to the VM

var_list = ['pr','tas']         #['tas','pr']
season_list = ['wet','dry'] 
base_period = [1991,2020]

# Parameters to change depending on region and time horizon
region_list = ['BES']  
target_period_ec_G = [2026,2055] # [2036,2065]
target_period_ec_W = [2046,2075] # [2081,2110]
target_period_cmip6 = [2035,2065] # [2070,2100]  
target_period_title = '2050' # '2085' 
timeseriesdir = '/data/brotons/resampling/Method-2014/esmvaltool_output/recipe_BES_2050_20220316_151634/preproc/local_resampling/' #'/data/brotons/resampling/Method-2014/esmvaltool_output/recipe_BES_2085_20220316_090941/preproc/local_resampling/' 

change = 'relative'        # relative or absolute for ev or pr, tas is always absolute
set_y_lim = False           # fix y-axis range if True
y_lim_lower = -1
y_lim_upper = 0.1

import getpass
import glob
import os
import subprocess
import sys

user=getpass.getuser()
figdir = '/data/brotons/resampling/Method-2014/output_plots/'

##---------- FUNCTIONS ----------##

import xarray as xr
import numpy as np
from datetime import datetime
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

tas_carib =  xr.open_mfdataset('/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/ECE_bis/tas_BES/tas_BES_region_mean/tas_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble',combine = 'nested')
pr_carib =  xr.open_mfdataset('/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/ECE_bis/pr_BES/pr_BES_region_mean/pr_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble',combine = 'nested')

ec_earth = {'pr':pr_carib , 'tas':tas_carib}


def time_period(da,period):
    """
    Select only the months inside the season of interest ('period')
    """
    if period == 'annual':
        da_time = da
    elif period == 'wet':
        da_time = da.where((da.time.dt.month>=5) & (da.time.dt.month<=11))
    elif period == 'dry':
        da_time = da.where((da.time.dt.month<=4) | (da.time.dt.month>=12))
    return da_time.dropna(dim='time')

def compute_change_MQ(file_in, season):
    """
    Open data, select data of interest (season, base&target period), 
    compute change (mean & quantiles),
    """
    # open data
    ds_in = xr.open_dataset(file_in)
    da_ts = ds_in[var]
    ds_in.close()

    # select months in season of interest
    da_ts = time_period(da_ts,season)
    da_ts.attrs['season'] = 'only months in '+season

    # split two periods
    da_base = da_ts.sel(time=slice(str(base_period[0])+'-01-01',str(base_period[1])+'-12-30'))
    da_target = da_ts.sel(time=slice(str(target_period_cmip6[0])+'-01-01',str(target_period_cmip6[1])+'-12-30'))
    da_base.attrs['period'] = '01-01-'+str(base_period[0])+' to 31-12-'+str(base_period[1])
    da_target.attrs['period'] = '01-01-'+str(target_period_cmip6[0])+' to 31-12-'+str(target_period_cmip6[1])

    # compute changes
    da_base_mean = da_base.mean()
    da_target_mean = da_target.mean()
    da_base_quant = da_base.quantile([0.05,0.1,0.5,0.9,0.95])
    da_target_quant = da_target.quantile([0.05,0.1,0.5,0.9,0.95])
    if change == 'absolute' or var == 'tas':
        da_change_mean = da_target_mean - da_base_mean
        da_change_quant = da_target_quant - da_base_quant
        if var in ['ev','pr']:
            if da_ts.attrs['units'] == 'kg m-2 s-1':
                da_change_mean = da_change_mean * 86400
                da_change_mean.attrs['units'] = 'mm day-1'  # mm/day
                da_change_quant = da_change_quant * 86400
                da_change_quant.attrs['units'] = 'mm day-1'  # mm/day
            else:
                print("Expected kg m2 s-1 as units for ev or pr")
                exit(1)
    elif change == 'relative' and var != 'tas':
        da_change_mean = (da_target_mean - da_base_mean)/da_base_mean*100
        da_change_quant = (da_target_quant - da_base_quant)/da_base_quant*100
        da_change_mean.attrs['units'] = '%'
        da_change_quant.attrs['units'] = '%'
    else: 
        print("Set change to 'absolute' or 'relative'")
        exit(1)
    return(da_change_mean, da_change_quant)

def print_model_list(filename, var):
    """
    Create a documentation file with all model-data taken into account in the figure
    """
    f = open(filename, "w")
    print(f"File created at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", file=f)
    print(f"With script band_resamples.py with git-revision {subprocess.check_output(['git', '--no-pager', 'show', '--date=short', '--format=%ad', '--name-only']).split()[0].decode() + '_' + subprocess.check_output(['git', '--no-pager', 'describe', '--tags', '--always', '--long', '--dirty']).strip().decode()}", file=f)
    print('CMIP6 ensemble:', file=f)
    os.listdir(f"{timeseriesdir}{var}_cmip/")
    f.close()

##---------- COMPUTE ----------##

fig,ax = plt.subplots(2,2,figsize=(15,10)) #to plot; if this is inside the loop only one subplot appears

for var in var_list:
    for season in season_list:
        # for ssp in ssp_list:
            for region in region_list:
                if region not in ['BES','BES Bonaire', 'BES St. Eustatius and Saba']:
                    print('Choose region BES,BES Bonaire, BES St. Eustatius and Saba')
                    exit(1)
                if var not in ['tas','pr']:
                    print('Choose var [tas,pr]')
                    exit(1)
                # if ssp not in ['126','245','585']:
                #     print('Choose ssp [126,245,585]')
                #     exit(1)
                if season not in ['dry','wet','annual']:
                    print('Choose season [dry,wet,annual]')
                    exit(1)

                datadir = timeseriesdir + region + '/' + var + '/'

                # one run per CMIP6 model
                model_mean_data = []
                model_quant_data = []
                for file_in in glob.iglob(f"{timeseriesdir}{var}_cmip/CMIP6*.nc"):
                    mean_data, quant_data = compute_change_MQ(file_in, season) 
                    model_mean_data.append(mean_data)
                    model_quant_data.append(quant_data.rename({'quantile': 'percentile'}))
                all_data_mean = xr.concat(model_mean_data, dim='model')
                all_data_quant = xr.concat(model_quant_data, dim='model')
                

                # compute multi model statistics
                mm_change_mean = all_data_mean.quantile((0.1,0.25,0.5,0.75,0.9),dim='model')
                mm_change_mean = mm_change_mean.rename({'quantile':'MMquant'})
                mm_change_quant = all_data_quant.quantile((0.1,0.25,0.5,0.75,0.9),dim='model')
                mm_change_quant = mm_change_quant.rename({'quantile':'MMquant'})

                # EC-Earth ensembles:

                ec_earth_data = ec_earth[var] #select variable

                #select season
                if season == 'dry': ec_earth_data = ec_earth_data.where((ec_earth_data.time.dt.month<=4) | (ec_earth_data.time.dt.month>=12))
                elif season == 'wet': ec_earth_data = ec_earth_data.where((ec_earth_data.time.dt.month>=5) & (ec_earth_data.time.dt.month<=11))

                #select scenario (W,G) depending on period
                ec_earth_G = ec_earth_data.sel(time=slice(str(target_period_ec_G[0]),str(target_period_ec_G[1])))
                ec_earth_W = ec_earth_data.sel(time=slice(str(target_period_ec_W[0]),str(target_period_ec_W[1])))
                ec_earth_base = ec_earth_data.sel(time=slice(str(base_period[0])+'-01-01',str(base_period[1])+'-12-30'))

                #compute statistics
                ec_earth_W_quant = ec_earth_W.quantile([0.05,0.1,0.5,0.9,0.95],dim='time')
                ec_earth_W_mean = ec_earth_W.mean(dim='time')
                ec_earth_G_quant = ec_earth_G.quantile([0.05,0.1,0.5,0.9,0.95],dim='time')
                ec_earth_G_mean = ec_earth_G.mean(dim='time')
                ec_earth_base_quant = ec_earth_base.quantile([0.05,0.1,0.5,0.9,0.95],dim='time')
                ec_earth_base_mean = ec_earth_base.mean(dim='time')

                if var == 'pr': 
                    ec_earth_W_change_quant = ((ec_earth_W_quant.pr - ec_earth_base_quant.pr) /ec_earth_base_quant.pr)*100
                    ec_earth_G_change_quant = ((ec_earth_G_quant.pr - ec_earth_base_quant.pr) /ec_earth_base_quant.pr)*100
                    ec_earth_W_change_mean = ((ec_earth_W_mean.pr - ec_earth_base_mean.pr) /ec_earth_base_mean.pr)*100
                    ec_earth_G_change_mean = ((ec_earth_G_mean.pr - ec_earth_base_mean.pr) /ec_earth_base_mean.pr)*100
                elif var=='tas': 
                    ec_earth_W_change_quant = ec_earth_W_quant.tas - ec_earth_base_quant.tas
                    ec_earth_G_change_quant = ec_earth_G_quant.tas- ec_earth_base_quant.tas
                    ec_earth_W_change_mean = ec_earth_W_mean.tas - ec_earth_base_mean.tas
                    ec_earth_G_change_mean = ec_earth_G_mean.tas - ec_earth_base_mean.tas

                ##---------- PLOT ----------##

                if var == 'pr': row = 0
                elif var == 'tas': row = 1
                if season == 'dry': column = 0; seasons = 'Dry season ';months = '(December-April)'
                elif season == 'wet': column = 1; seasons = 'Wet season ';  months = '(May-November)'

                #plt.subplots_adjust(left=0.14, bottom=0.12, right=.95, top=0.92)
                # CMIP6 band: mean
                ax[row,column].fill_between([-.5,.5],mm_change_mean.sel(MMquant=0.1),mm_change_mean.sel(MMquant=0.9),color='lightgrey',linewidth=0)
                ax[row,column].fill_between([-.5,.5],mm_change_mean.sel(MMquant=0.25),mm_change_mean.sel(MMquant=0.75),color='grey',linewidth=0)
                ax[row,column].plot([-.5,.5],np.repeat(mm_change_mean.sel(MMquant=0.5).values,2),color='black',linewidth=1.5)
                # CMIP6 band: quant
                ax[row,column].fill_between(range(1,6),mm_change_quant.sel(MMquant=0.1),mm_change_quant.sel(MMquant=0.9),color='lightgrey',linewidth=0)
                ax[row,column].fill_between(range(1,6),mm_change_quant.sel(MMquant=0.25),mm_change_quant.sel(MMquant=0.75),color='grey',linewidth=0)
                ax[row,column].plot(range(1,6),mm_change_quant.sel(MMquant=0.5),color='black',linewidth=1.5)

                for i in range(0,len(tas_carib.ensemble)):

                #EC-Earth : 
                    ax[row,column].plot(range(1,6),ec_earth_W_change_quant[:,i,0,0],color='orange',linewidth=1.5)
                    ax[row,column].plot(range(1,6),ec_earth_G_change_quant[:,i,0,0],color='blue',linewidth=1.5)
                    ax[row,column].plot([-.5,.5],[ec_earth_W_change_mean[i,0,0],ec_earth_W_change_mean[i,0,0]],color='orange',linewidth=1.5)
                    ax[row,column].plot([-.5,.5],[ec_earth_G_change_mean[i,0,0],ec_earth_G_change_mean[i,0,0]],color='blue',linewidth=1.5)

                # # labels etc
                plt.setp(ax, xticks=range(0,6), xticklabels=['ave','P05','P10','P50','P90','P95'])
                if var == 'tas':
                     ax[row,column].set_title(seasons+months,loc='right')
                     ax[row,column].set_ylabel(r'$\Delta$ Temperature [$^\circ$ C]',fontsize=14)
                elif var == 'pr' and change == 'relative':
                     ax[row,column].set_title(seasons+months,loc='right')
                     ax[row,column].set_ylabel(r'$\Delta$ Precipitation [%]',fontsize=14)
                elif var == 'pr' and change == 'absolute':
                     ax[row,column].set_title(seasons+months,loc='right')
                     ax[row,column].set_ylabel(r'$\Delta$ Precipitation [mm/day]',fontsize=14)
                
                plt.suptitle(f"{region_list[0]} region; {target_period_title}",fontsize=18)

                cmip6_80 =  mpatches.Patch( color='lightgrey', label=r'$80\%$ CMIP6 ')
                cmip6_50 = mpatches.Patch(color='grey', label=r'$50\%$ CMIP6 ')
                ecearth_W = mlines.Line2D([], [], color='orange', label=r'W ($dT = 1 ^{\circ} C$)')
                ecearth_G = mlines.Line2D([], [], color='blue', label=r'G ($dT = 2 ^{\circ} C$)')
                plt.legend(handles=[cmip6_80,cmip6_50,ecearth_W,ecearth_G])
                
                output_filename = figdir+region+'_'+target_period_title +'_fig4_Lenderink2014.png'
                fig.savefig(output_filename)





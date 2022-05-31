# Python script to differentate between dry and wet models. Output: a figure and the list of wet and dry models. (Figure 21 report)
# Wet models are the models with a scaled change in precipitation averaged over the Caribbean region wetter than the 70th percentile
# Dry models are the models with a scaled change in precipitation averaged over the Caribbean region dryier than the 30 th percentile


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import os
import pandas as pd
import datetime
from cftime import DatetimeNoLeap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
from matplotlib import cm

# ------------------ Things to select ----------------------

scaling = True
ssp = 'ssp585'
season_list = ['annual']#,'dry','wet']

# ----------------------------------------------------------


#Functions

def time_period(da,period):
    """
    Select only the months inside the season of interest ('period') and takes running mean 30 yr window
    """
    if period == 'annual':
        da_time = da
        #da_time = da_time.rolling(time=30*12,center=False).mean(dim='time')
    elif period == 'wet':
        da_time = da.where((da.time.dt.month>=5) & (da.time.dt.month<=11),drop=True)#.dropna(dim='time')
        #da_time = da_time.rolling(time=30*7,center=False).mean(dim='time')
    elif period == 'dry':
        da_time = da.where((da.time.dt.month<=4) | (da.time.dt.month>=12),drop=True)#.dropna(dim='time')
        #da_time = da_time.rolling(time=30*5,center=False).mean(dim='time')
    return da_time #.dropna(dim='time')

def boxmean(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.lat))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('lat','lon'))
    return boxmean


def preprocessing_area_mean(ds):
    """
    Function to call when importing netcdf files with xr.open_mfdataset.
    - Sets all the datasets to the same datetime and time axis
    - Deletes from the files inconsistent dimensions/coordinates and saves models name
    """ 
    if ds.indexes['time'].inferred_type != "datetime64":
            ds['time'] = ds.indexes['time'].to_datetimeindex()
    ds['time'] = pd.date_range(start=str(ds.time[0].dt.year.values)+'-0'+str(ds.time[0].dt.month.values)+'-01',end='2100-12-16',periods=len(ds.time))
    if 'plev_bnds' in ds:
        ds =  ds.drop('plev_bnds')
    if 'plev' in ds:
        ds = ds.drop('plev')
    if 'height' in ds:
        ds = ds.drop('height')
    if 'file_qf' in ds:
        ds = ds.drop('file_qf')
    return ds

#------------------------INPUT------------------------#

season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}

#------------------------CALCULATIONS------------------------#

# Import data

pr_BES_cmip6 = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES2/pr/*'+ssp+'*.nc',concat_dim = 'model',combine = 'nested',preprocess = preprocessing_area_mean)
tas_40S_40N_cmip6 = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*'+ssp+'*.nc',concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
pr_BES_cmip6 = pr_BES_cmip6.load()
tas_40S_40N_cmip6 = tas_40S_40N_cmip6.load()

pr_BES_cmip6 = pr_BES_cmip6.sortby(pr_BES_cmip6.model_name) # to avoid different models order!
tas_40S_40N_cmip6 = tas_40S_40N_cmip6.sortby(tas_40S_40N_cmip6.model_name)

models = tas_40S_40N_cmip6.model_name




#---- Define wet and dry models:-----

for season in season_list:

    tas_period = time_period(tas_40S_40N_cmip6,season)
    pr_period = time_period(pr_BES_cmip6,season)

    # Compute change in tempertature
    dtas_40S_40N = tas_period.sel(time=slice('1-1-2071', '31-12-2100')).mean(dim='time') - tas_40S_40N_cmip6.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') # target-ref
    dpr_area_mean = pr_period.sel(time=slice('1-1-2071', '31-12-2100')).mean(dim='time') - pr_BES_cmip6.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') # target-ref
    
    if scaling == True:
        dpr_dt40 = dpr_area_mean.pr/dtas_40S_40N.tas #scaled region mean precipitation
        ylaebl_name = r'Scaled $\Delta$ Precipitation (Carib) [mm/day $^{\circ}$C]'; sc = ''; title = 'Scaled precipitation'
    elif scaling == False:
        dpr_dt40 = dpr_area_mean.pr
        ylaebl_name = r'$\Delta$ Precipitation (Carib) [mm/day $^{\circ}$C]'; sc = '_nonscaled'; title = 'Non-scaled precipitation'
    dpr_dt40_quant = dpr_dt40.quantile((0.3,0.7))
    # dpr_dt40_quant = dpr_dt40.quantile((0.45,0.55))
    wet_idx = np.where(dpr_dt40>=dpr_dt40_quant[1])[0]
    dry_idx = np.where(dpr_dt40<=dpr_dt40_quant[0])[0]

    print('Scaling = '+str(scaling))
    print(season)
    models_wet = []
    print('Wet models are:')
    for i in range(len(wet_idx)):
        models_wet.append(models[int(wet_idx[i])])
        print('-', models[int(wet_idx[i])].values)


    models_dry = []
    print('Dry models are:')
    for i in range(len(dry_idx)):
        models_dry.append(models[int(dry_idx[i])])
        print('-', tas_40S_40N_cmip6.model_name[int(dry_idx[i])].values) 

    models_names = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'AWI-CM-1-1-MR', 'BCC-CSM2-MR','    CESM2', 'CESM2-WACCM   ', '    CIESM', 'CMCC-CM2-SR5', 'CNRM-CM6-1','CNRM-CM6-1-HR', 'CNRM-ESM2-1', 'CanESM5', 'CanESM5-CanOE','EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'FGOALS-g3', 'GFDL-ESM4', 'GISS-E2-1-G', 'HadGEM3-GC31-LL   ', '     INM-CM4-8   ', '   INM-CM5-0  ', '   IPSL-CM6A-LR    ', 'KACE-1-0-G', 'MIROC-ES2L', 'MIROC6','MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-MM', 'UKESM1-0-LL']
    breakpoint()
    plt.figure(345678)
    fig,ax = plt.subplots(1,1,figsize=(15,10))
    #fig.suptitle(title)
    ax.set_title(season_title[season],loc='right',fontsize=16)
    ax.set_title(title_ssp[ssp],loc='left',fontsize=16)
    plt.scatter(range(0,len(models)),dpr_dt40* 86400,s=100,c='black')
    plt.ylabel(ylaebl_name,fontsize=18)
    plt.axhline(dpr_dt40_quant[0]* 86400, color='peru', linestyle='--',label='P30')
    plt.axhline(dpr_dt40_quant[1]* 86400, color='teal', linestyle='--',label='P70')
    plt.scatter(dry_idx,dpr_dt40[dry_idx]* 86400,s = 200, c = 'peru')
    plt.scatter(wet_idx,dpr_dt40[wet_idx]* 86400,s = 200, c = 'teal')
    ax.set_xticks(np.linspace(-1.5,30.5,len(models)))#,rotation=50)
    ax.set_xticklabels(models_names,rotation=50,fontsize = 14)
    #ax.set_xticklabels(np.arange(1,33),rotation=50,fontsize = 14)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16)
    plt.tight_layout()
    # plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/wet_dry_deffinition_'+season+sc+'.png')
    plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/wet_dry_deffinition_'+season+sc+'.png',dpi = 1200)
# This script computes correlation maps. Correlation between scaled change in temperature (or wap) and scaled change in precipitation averaged over the Caribbean region (Figures 20 and A6)

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import os
import pandas as pd
import datetime
from cftime import DatetimeNoLeap
from tqdm import tqdm
#from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import cm
from scipy import stats

import statsmodels.api as sm
import numpy as np
from scipy.stats import t

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#------------------------INPUT------------------------#

# SETTINGS PROGRAM

# -------------Things to choose:----------------------

var = 'tas' #wap500 or psl
pr_units = 'mm/day' #mm/day or %
ssp = 'ssp585'
confidence_level = 0.9

domain = 'NAtl_NPac'

# ----------------------------------------------------

# list things to do
season_list = ['dry','wet']
which_models_list = ['all']#,'dry','wet']

# Paths variables
paths = {'ua850':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/ua850/', 
        'va850':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/va850/',
        'ua200':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/ua200/' ,
        'va200':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/va200/',
        'pr': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/pr/',
        'hfls': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/hfls/',
        'hus850': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/hus850/',
        'psl': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/psl/',
        'tas': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/tas/',
        'wap500': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/wap500/',
        'moist_conv': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_conv/',
        'moist_mass_conv':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_mass_conv/',
        'moist_adv': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_adv/'}

path_area_mean = {'tas_40': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*'+ssp+'*.nc',
                'pr_BES':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/pr/*'+ssp+'*.nc' }

# Wet and dry models for scaled and non scaled precipitation
wet_models =['CIESM', 'CNRM-CM6-1-HR', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'FGOALS-f3-L', 'INM-CM4-8', 'INM-CM5-0', 'NESM3'] #for annnual and scaled precipitation!
dry_models = ['ACCESS-ESM1-5', 'CESM2', 'CanESM5-CanOE', 'CanESM5', 'FGOALS-g3', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR']
dry_wet_models = {'wet':wet_models,'dry':dry_models}

# definition coordinates
latit = {'BES': [5,20], 'NAtl_NPac': [-20,70],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40],'Atl_Pac': [-90,90]}
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[120,170],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340],'Atl_Pac': [110,360]}

latitude = latit[domain]
longitude = longit[domain]

# Pre-settings (they can be changed along the code):
gr = ''

# SETTINGS FIGURES

season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}
title = {'tas':r'Correlation BES sc. $\Delta$ precipitation and sc. $\Delta$ temperature change (2005-2075). ','psl':r'Correlation BES sc. $\Delta$ precipitation and sc. $\Delta$ sea level pressure change (2005-2075). ','wap500':r'Correlation BES sc. $\Delta$ precipitation and sc. $\Delta$ Omega at 500hPa change (2005-2075).' }
title_confidence = str(confidence_level) + 'confidence range. '
# title = {'tas':r'Slope BES sc. $\Delta$ precipitation and sc. $\Delta$ temperature change (2005-2075). ','psl':r'Slope BES sc. $\Delta$ precipitation and sc. $\Delta$ sea level pressure change (2005-2075). ','wap500':r'Slope BES sc. $\Delta$ precipitation and sc. $\Delta$ Omega at 500hPa change (2005-2075).' }
title_which_models = {'all': '','dry': 'Dry models. ', 'wet': 'Wet models. ', 'dry-wet': 'Dry-Wet models. '}
cbar_title = {'tas': r'Corr. coef. (Sc $\Delta$ Temp and Sc $\Delta$ Pr BES)','psl':'Corr. coef. (Sc $\Delta$SLP and Sc $\Delta$ Pr BES)','wap500':'Corr. coef. (Sc $\Delta$ Omega and Sc $\Delta$ Pr BES)'}


#Functions

def time_period(da,period):
    """
    Select only the months inside the season of interest ('period') and takes running mean 30 yr window
    """
    if period == 'annual':
        da_time = da
    elif period == 'wet':
        da_time = da.where((da.time.dt.month>=5) & (da.time.dt.month<=11),drop=True)
    elif period == 'dry':
        da_time = da.where((da.time.dt.month<=4) | (da.time.dt.month>=12),drop=True)
    return da_time 

def time_month(da,month, n_months):
    """
    Select only the months inside the season of interest ('month') and takes running mean 30 yr window
    """

    if n_months == 1:
        da_time = da.where(da.time.dt == month, drop=True)
    elif n_months >= 2:  
        if month[0] < month[-1]:
            da_time = da.where((da.time.dt.month>=month[0]) & (da.time.dt.month<=month[-1]),drop=True)
        elif month[0]>month[-1]:
            da_time = da.where((da.time.dt.month==month[0]) | (da.time.dt.month==month[-1]),drop=True)
    return da_time 

def boxmean(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.lat))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('lat','lon'))
    return boxmean

def preprocessing_lat_lon(ds):
    """
    Function to call when importing netcdf files with xr.open_mfdataset.
    - Sets all the datasets to the same datetime and time axis
    - Selects the domain we are interested on
    - Deletes from the files inconsistent dimensions/coordinates and saves models name
    - Saves the name of the model and the rpf number
    """
    ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]) & (ds.lon > longitude[0]) & (ds.lon < longitude[1]), drop=True)   
    if 'plev_bnds' in ds:
        ds =  ds.drop('plev_bnds')
    if 'plev' in ds:
        ds = ds.drop('plev')
    if 'height' in ds:
        ds = ds.drop('height')
    if 'file_qf' in ds:
        ds = ds.drop('file_qf')
    return ds

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


def test_cor(x,y,conf):
    """
    Function provided by Simon Michel (IMAU).
    To compute significance of the correlation coefficient with the mehodology from McCarthy et al. 2015
    """

    ax=sm.tsa.acf(x)[1]
    ay=sm.tsa.acf(y)[1]
   
    r=np.corrcoef(x,y)[0,1]
    N=np.size(x)

    Neff=N*((1-ax*ay)/(1+ax*ay))

    #print(r)
    #print(Neff)
   
    TT=np.sqrt(Neff)*(r/np.sqrt(1-r**2))

    if abs(TT)>t.interval(conf+(1-conf)/2, Neff, loc=0, scale=1)[1]:
        return(True)
    else:
        return(False)



#------------------------CALCULATIONS------------------------#
breakpoint()
# Import variables
# precipitation
pr_BES = xr.open_mfdataset(path_area_mean['pr_BES'],concat_dim = 'model',combine = 'nested',preprocess = preprocessing_area_mean)
pr_BES = pr_BES.load()
pr_BES = pr_BES.pr * 86400
pr_BES = pr_BES.sortby(pr_BES.model_name)


# temperature to scale
tas_40 = xr.open_mfdataset(path_area_mean['tas_40'],concat_dim = 'model',combine = 'nested',preprocess = preprocessing_area_mean)
tas_40 = tas_40.load()
tas_40 = tas_40.tas
tas_40 = tas_40.sortby(tas_40.model_name)

# Calculations

models = tas_40.model_name

for season in season_list:

    if var == 'tas':
        tas_2005 = xr.open_mfdataset(paths['tas']+'*_'+season+'_season_2005.nc',concat_dim = 'model',combine = 'nested', preprocess = preprocessing_lat_lon)
        tas_2085 = xr.open_mfdataset(paths['tas']+'*_'+season+'_season_2085.nc',concat_dim = 'model',combine = 'nested', preprocess = preprocessing_lat_lon)
        tas_2005 = tas_2005.load(); tas_2085 = tas_2085.load()
        ds_2005 = tas_2005.tas
        ds_2085 = tas_2085.tas

    elif var == 'wap500':
        wap500_2005 = xr.open_mfdataset(paths['wap500']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        wap500_2085 = xr.open_mfdataset(paths['wap500']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        wap500_2005 = wap500_2005.load(); wap500_2085 = wap500_2085.load()
        ds_2005 = wap500_2005.wap[:,0,:,:]
        ds_2085 = wap500_2085.wap[:,0,:,:]

    elif var1 == 'psl':
        psl_2005 = xr.open_mfdataset(paths['psl']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        psl_2085 = xr.open_mfdataset(paths['psl']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        psl_2005 = psl_2005.load(); psl_2085 = psl_2085.load()
        ds_2005 = psl_2005.wap[:,0,:,:]
        ds_2085 = psl_2085.wap[:,0,:,:]

    ds_2005 = ds_2005.sortby(ds_2005.model_name)
    ds_2085 = ds_2085.sortby(ds_2085.model_name)

    pr_period = time_period(pr_BES,season)
    #tas_40_period = time_period(tas_40_models,season)

    #Caribbean precipitation
    pr_2005 = pr_period.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time')
    pr_2085 = pr_period.sel(time=slice('1-1-2071', '31-12-2100')).mean(dim='time')

    #Temperature
    # tas_40_2005 = tas_40_period.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time')
    # tas_40_2085 = tas_40_period.sel(time=slice('1-1-2071', '31-12-2100')).mean(dim='time')

    tas_40_2005 = tas_40.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time')
    tas_40_2085 = tas_40.sel(time=slice('1-1-2071', '31-12-2100')).mean(dim='time')

    
    ds = (ds_2085-ds_2005)/(tas_40_2085-tas_40_2005)
    pr = (pr_2085-pr_2005)/(tas_40_2085-tas_40_2005)

    for which_models in which_models_list:

        if which_models == 'all':
            ds_models = ds
            pr_models = pr

        elif which_models == 'dry' or which_models == 'wet':
            models = dry_wet_models[which_models]
            models_idx = np.zeros(len(models))
            for i in range(0,len(models)):
                models_idx[i] = np.where(ds.model_name == models[i])[0]

            models_idx = models_idx.astype(int)

            ds_models = ds.sel(model = models_idx)
            pr_models = pr.sel(model = models_idx)



        # compute correlation coefficient bewteen pr_BES and every gridpoint od the second var
        corrcoef = np.zeros((len(ds_models.lat),len(ds_models.lon)))
        confidence = np.zeros((len(ds_models.lat),len(ds_models.lon)))
        
        breakpoint()
        for i in range(0,len(ds_models.lat)):
            for j in range(0,len(ds_models.lon)):
                idx = np.isfinite(pr_models[:]) & np.isfinite(ds_models[:,i,j])
                corrcoef[i,j] = (np.corrcoef(pr_models[idx],ds_models[idx,i,j])[0,1])

                
                if test_cor(pr_models[idx],ds_models[idx,i,j],confidence_level):
                    confidence[i,j] = 1
    
        breakpoint()
        # to set corrcoef same coordinates as tas
        ds_coords = ds_models[0,:,:]
        corrcoef_toplot = xr.DataArray(corrcoef, coords=ds_coords.coords)
        confidence_toplot = xr.DataArray(confidence, coords=ds_coords.coords)


        #plot figure 
        plt.rcParams['hatch.linewidth'] = 0.05
        fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
        fig.suptitle(title[var]+title_which_models[which_models]+title_confidence)
        ax.set_title(title_ssp[ssp],loc='left',fontsize=14) #12
        ax.set_title(season_title[season],loc='right',fontsize=14) #14
        z1_plot = ax.pcolormesh(corrcoef_toplot.lon,corrcoef_toplot.lat,corrcoef_toplot,transform=ccrs.PlateCarree(central_longitude=0), cmap=matplotlib.cm.get_cmap('cmo.curl') ) 
        zconf_plot = plt.contourf(confidence_toplot.lon, confidence_toplot.lat, confidence_toplot, levels=[0.99, 2],colors='none', hatches=['...'], transform=ccrs.PlateCarree(central_longitude= 0))
        ax.coastlines(resolution = '50m')
        ax.set_xlim([-70,180])
        ax.set_ylim([-20,70]) #40
        ax.plot([90,90],[5,20],color = 'red',linewidth='1')
        ax.plot([130,130],[5,20],color = 'red',linewidth='1')
        ax.plot([90,130],[5,5],color = 'red',linewidth='1')
        ax.plot([90,130],[20,20],color = 'red',linewidth='1')

        # ax.plot([90,90],[5,20],color = 'black',linewidth='1')
        # ax.plot([130,130],[5,20],color = 'black',linewidth='1')
        # ax.plot([90,130],[5,5],color = 'black',linewidth='1')
        # ax.plot([90,130],[20,20],color = 'black',linewidth='1')

        gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude= 0), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 11  , 'color': 'black'}
        gl.ylabel_style = {'size': 11, 'color': 'black'}

        axins = inset_axes(ax,width="70%",  height="8%",loc='lower center',borderpad=-4)

        z1_plot.set_clim((-1,1))
        cbar = fig.colorbar(z1_plot, cax=axins, shrink=0.5, orientation='horizontal')
        cbar.set_label(cbar_title[var],fontsize = 14)
        cbar.ax.tick_params(labelsize = 11)

        # fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/correlation_map_pr_'+var+gr+'_'+season+'_season_'+which_models+'_models_'+ssp+'.png')
        #fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/pwpt/correlation_map_pr_'+var+gr+'_'+season+'_season_'+which_models+'_models_'+ssp+'.png')
        fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/correlation_map_pr_'+var+gr+'_'+season+'_season_'+which_models+'_models_'+ssp+'.png',dpi = 300)
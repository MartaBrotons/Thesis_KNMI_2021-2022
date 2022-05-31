# This scripts plots for precipitation and temperature at 2 m whisker plots for the 2036-2065 and 2071-2100 averages. (Figure A.1)

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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import seaborn as sns


#------------------------INPUT------------------------#

season_list = ['annual']#,'dry','wet']
ssp_list = ['ssp585']#,'ssp245','ssp126']
var_list = ['tas','pr']

season_title = {'annual':'Annual','dry':'Dry season (December-April)','wet':'Wet season (May-November)'}
ssp_title = {'ssp585':'SSP5-8.5','ssp245':'SSP2-4.5','ssp126':'SSP1-2.6'}
colors = {'tas':'r','pr':'b'}
y_label = {'tas':r'$\Delta$ Temp. [$^{\circ}C$/$^{\circ}C$] (w.r.t. 1991-2020)','pr': r'$\Delta$ Precip. [%/$^{\circ}C$] (w.r.t. 1991-2020)'}

#color = {'tas':matplotlib.cm.get_cmap('cmo.curl'), 'psl':cm.PuOr}
pr_units = '%' #mm/day or %


paths_pr = {'ssp585':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/pr/*ssp585*.nc',
            'ssp245':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/pr/*ssp245*.nc'}
paths_tas = {'ssp585':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/tas/*ssp585*.nc',
            'ssp245': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/tas/*ssp245*.nc'}

paths_tas_40 = {'ssp585':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*ssp585*.nc',
            'ssp245': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*ssp245*.nc'}

path = {'tas': paths_tas, 'pr': paths_pr}

#------------------------CALCULATIONS------------------------#

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

def preprocessing(ds):
    'Function to call when importing netcdf files with xr.open_mfdataset.'
    'Deletes from the files inconsistent dimensions/coordinates and saves models name'
    if 'plev_bnds' in ds:
        ds =  ds.drop('plev_bnds')
    if 'plev' in ds:
        ds = ds.drop('plev')
    if 'height' in ds:
        ds = ds.drop('height')
    if 'file_qf' in ds:
        ds = ds.drop('file_qf')
    return ds.assign_coords(model_name = ds.source_id)

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
# --------------------------------------------------------------------------#

breakpoint()
for season in tqdm(season_list):
    for ssp in ssp_list:

        tas_BES_40 = xr.open_mfdataset(paths_tas_40[ssp],concat_dim = 'model',combine = 'nested',preprocess = preprocessing_area_mean) #import data
        tas_BES_40 = tas_BES_40.load()

        for var in var_list:
        
            # Import data
            ds_BES = xr.open_mfdataset(path[var][ssp],concat_dim = 'model',combine = 'nested',preprocess = preprocessing_area_mean) #import data
            ds_BES = ds_BES.load()

            ds_period = time_period(ds_BES,season) #period selection annual, dry, wet
            tas_BES_40_period = time_period(tas_BES_40,season)

            ds_ref = ds_period.sel(time=slice('01-01-1991', '31-12-2020')).mean(dim='time') #reference mean
            tas_BES_40_ref = tas_BES_40.sel(time=slice('01-01-1991', '31-12-2020')).mean(dim='time')

            #change in precipitation
            if var == 'pr':
                if pr_units == '%':
                    dds = (ds_period-ds_ref)*100/ds_ref
                    y_label['pr'] = r'$\Delta$ Precip. [%/ $^{\circ}C$] (w.r.t. 1991-2020)'
                    pr_name = '%'
                elif pr_units =='mm/day':
                    dds = (ds_period-ds_ref) * 86400
                    y_label['pr'] = r'$\Delta$ Precip. [mm/day $^{\circ}C$] (w.r.t. 1991-2020)'
                    pr_name = 'mmday'
                dds_scaled = dds.pr#/(tas_BES_40_period.tas-tas_BES_40_ref.tas)
                dds_MOC = dds_scaled.sel(time=slice('01-01-2035', '31-12-2065')).mean(dim='time')#.pr
                dds_EOC = dds_scaled.sel(time=slice('01-01-2070', '31-12-2100')).mean(dim='time')#.pr
            elif var == 'tas': 
                dds_scaled = (ds_period.tas-ds_ref.tas)#/(tas_BES_40_period.tas-tas_BES_40_ref.tas)
                dds_MOC = dds_scaled.sel(time=slice('01-01-2035', '31-12-2065')).mean(dim='time')#.tas
                dds_EOC = dds_scaled.sel(time=slice('01-01-2070', '31-12-2100')).mean(dim='time')#.tas

            dds_boxplot = [dds_MOC,dds_EOC]

            #plot figure 
            fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True)
            ax.set_title(season_title[season],loc='right',fontsize=14)
            ax.set_title(ssp_title[ssp],loc='left',fontsize=14)
            ax.set_ylabel(y_label[var],fontsize=14)
            #ax.set_xlabel('Time [year]',fontsize=14)
            bplot1 = ax.boxplot(dds_boxplot,sym = '',showfliers=None, vert=True, patch_artist=True,boxprops=dict(facecolor=colors[var])) 
            ax.set_xticklabels(['2050 (2035-2065)','2085 (2070-2100)'],fontsize=14)
            #sns.boxplot(x='variable', y='value', data=dpr_boxplot)

            #fig.savefig('/data/brotons/resampling/Method-2014/output_plots/correlation_map_pr_'+var+'_'+season+'_'+ssp+'.png')
            fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014//output_plots/report/CMIP6_boxplot_'+var+'_'+season+'_season_'+ssp+'.png',dpi = 1200)
            breakpoint()

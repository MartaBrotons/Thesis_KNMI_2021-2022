# This script plots precipitation change projections averaged over the Caribbean region for 32 models. Plots the multimode mean + 90% spread. (Figure 1)
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


#------------------------INPUT------------------------#

season_list = ['annual']#,'dry','wet']
ssp_list = ['ssp585']#,'ssp245','ssp126']
n_month = {'annual':12,'dry':5,'wet':7}
season_title = {'annual':'Annual','dry':'Dry season (December-April)','wet':'Wet season (May-November)'}
label = {'ssp585':'SSP5-8.5','ssp245':'SSP2-4.5','ssp126':'SSP1-2.6'}
colors = {'ssp585':'darkred','ssp245':'orange','ssp126':'midnightblue'}
colors_shaded = {'ssp585':'lightcoral','ssp245':'navajowhite','ssp126':'lightsteelblue'}

#color = {'tas':matplotlib.cm.get_cmap('cmo.curl'), 'psl':cm.PuOr}
color = {'tas':matplotlib.cm.get_cmap('cmo.matter'), 'psl':matplotlib.cm.get_cmap('cmo.turbid')}
pr_units = '%' #mm/day or %


# Working directories var: 
# wdir = {'ssp585':'/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/CMIP6/pr_BES/*ssp585*.nc',
#         'ssp245':'/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/CMIP6/pr_BES/*ssp245*.nc',
#         'ssp126':'/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/CMIP6/pr_BES/ssp126/*ssp126*.nc'}

wdir = {'ssp585':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/pr/*ssp585*.nc'}

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
    dpr_mean = []
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True)

    for ssp in ssp_list:

        pr_BES = xr.open_mfdataset(wdir[ssp],concat_dim = 'model',combine = 'nested', preprocess= preprocessing_area_mean) #import data
        pr_BES = pr_BES.load()
        pr_period = time_period(pr_BES,season) #period selection annual, dry, wet

        breakpoint()
        # pr_grouped = pr_period.groupby('time.year').mean(dim='time')
        # pr = pr_grouped.rolling(year=30,center=False).mean(dim='year').dropna('year') #rolling mean
        pr = pr_period.groupby('time.year').mean(dim='time')
        #pr = pr_grouped.rolling(year=30,center=False).mean(dim='year').dropna('year') #rolling mean
        pr_ref = pr.sel(year=slice('1991', '2020')).mean(dim='year') #reference mean

        breakpoint()
        #change in precipitation
        if pr_units == '%':
            dpr = (pr-pr_ref)*100/pr_ref
            Precip_label = r'$\Delta$ Precipitation [%]'
            pr_name = '%'
        elif pr_units =='mm/day':
            dpr = (pr-pr_ref) * 86400
            Precip_label = r'$\Delta$ Precipitation [mm/day]'
            pr_name = 'mmday'

        dpr = dpr.where(dpr.year>=1991)

        breakpoint()
        #plot figure 
        ax.set_title(season_title[season],loc='right',fontsize=13)
        ax.set_title('Precipitation Caribbean. SSP5-8.5 ',loc='left',fontsize=14)
        ax.set_ylabel(r'$\Delta$ Precip. [%] (w.r.t. 1991-2020)',fontsize=14)
        ax.axhline(0, color='black', linestyle='--')
        ax.set_xlabel('Year',fontsize=14)
        #ax.plot(dpr.time,dpr.pr.transpose(), color = 'grey')#colors_shaded[ssp]) 
        ax.fill_between(dpr.year,dpr.pr.quantile(0.05,dim='model'),dpr.pr.quantile(0.95,dim='model'),color='lightcoral',linewidth=0)
        #ax.plot(dpr.time,dpr.pr.mean(dim = 'model'), color = colors[ssp], label = label[ssp]) 
        ax.tick_params(axis='both', which='major', labelsize=14)
        dpr_mean.append(dpr.pr.mean(dim = 'model'))

    breakpoint()
    ax.plot(dpr.year,dpr_mean[0], color = colors['ssp585'], label = label[ssp],linewidth=2.0)  
    #ax.plot(dpr.time,dpr_mean[1], color = colors['ssp245'], label = label[ssp],linewidth=5.0)  
    #ax.plot(dpr.time,dpr_mean[1], color = colors['ssp126'], label = label[ssp],linewidth=5.0)  
    cmip6_90 =  mpatches.Patch( color='lightcoral', label=r'$90\%$ CMIP6 ')
    mean = mlines.Line2D([], [], color='darkred', label='CMIP6 multimodel mean')
    plt.legend(handles=[cmip6_90,mean],fontsize=13)

    #fig.savefig('/data/brotons/resampling/Method-2014/output_plots/correlation_map_pr_'+var+'_'+season+'_'+ssp+'.png')
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014//output_plots/report/CMIP6_pr_timeseries_'+season+'.png',dpi = 1200)

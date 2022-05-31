# This script computes the double ITCZ bias as CMIP6 - GPCP (1991-2020 averages). Then it computes the multimodel mean. (Figure 6)

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
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from scipy import stats
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#------------------------INPUT------------------------#


# SETTINGS PROGRAM

# -------------Things to choose:----------------------
# Variables

ssp = 'ssp585'
domain = 'NAtl_NPac'

# ----------------------------------------------------

# list things to do
season_list = ['annual']#,'dry','wet']

# defining domains
# definition coordinates
latit = {'BES': [5,20], 'NAtl_NPac': [-20,60],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40]}
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[120,170],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340]}

latitude = latit[domain]
longitude = longit[domain]

#for annual non-scaled precipitation
wet_models =['CIESM','CNRM-CM6-1-HR','CNRM-CM6-1','CNRM-ESM2-1','FGOALS-f3-L','INM-CM4-8','INM-CM5-0','NESM3'] 
dry_models =['ACCESS-ESM1-5','CESM2-WACCM','CESM2','CanESM5-CanOE','CanESM5','HadGEM3-GC31-LL','IPSL-CM6A-LR','KACE-1-0-G']
dry_wet_models = {'wet':wet_models,'dry':dry_models}

# Paths variables
paths = {'pr_cmip6': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/pr/*'+ssp+'*.nc',
        'pr_gpcp': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/gpcp_global_regridded.nc'}

# SETTINGS FIGURES

# Titles

title = 'CMIP6 bias: Multimodel mean - GPCP observations (1979-2020)'
season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}

# Functions

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

def preprocessing_lat_lon(ds):
    """
    Function to call when importing netcdf files with xr.open_mfdataset.
    - Sets all the datasets to the same datetime and time axis
    - Selects the domain we are interested on
    - Deletes from the files inconsistent dimensions/coordinates and saves models name
    - Saves the name of the model and the rpf number
    """
    if ds.indexes['time'].inferred_type != "datetime64":
            ds['time'] = ds.indexes['time'].to_datetimeindex()
    ds['time'] = pd.date_range(start=str(ds.time[0].dt.year.values)+'-0'+str(ds.time[0].dt.month.values)+'-01',end='2100-12-16',periods=len(ds.time))
    if domain == '40N_40S': ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]), drop=True)
    else: ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]) & (ds.lon > longitude[0]) & (ds.lon < longitude[1]), drop=True)   
    if 'plev_bnds' in ds:
        ds =  ds.drop('plev_bnds')
    if 'plev' in ds:
        ds = ds.drop('plev')
    if 'height' in ds:
        ds = ds.drop('height')
    if 'file_qf' in ds:
        ds = ds.drop('file_qf')
    ds = ds.assign_coords(model_name = ds.source_id)
    ds = ds.assign_coords(variant_label = ds.variant_label)
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

#------------------------CALCULATIONS------------------------#

# Import variables
breakpoint()

pr_cmip6 = xr.open_mfdataset(paths['pr_cmip6'],combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
pr_cmip6 = pr_cmip6.load()
pr_cmip6 = pr_cmip6.where((pr_cmip6.time.dt.year>=1991) & (pr_cmip6.time.dt.year<=2020),drop = True)
pr_cmip6 = pr_cmip6.pr*86400

pr_gpcp = xr.open_dataset(paths['pr_gpcp'])
pr_gpcp = pr_gpcp.where((pr_gpcp.lat > latitude[0]) & (pr_gpcp.lat < latitude[1]) & (pr_gpcp.lon > longitude[0]) & (pr_gpcp.lon < longitude[1]), drop=True)
pr_gpcp = pr_gpcp.where((pr_gpcp.time.dt.year>=1991) & (pr_gpcp.time.dt.year<=2020),drop=True) 
pr_gpcp = pr_gpcp.precip

#----------------------------------#

# Start computations

for season in season_list:
    breakpoint()

    pr_cmip6_season = time_period(pr_cmip6,season)
    pr_gpcp_season = time_period(pr_gpcp,season)

    pr_cmip6_mean = pr_cmip6_season.mean(dim = 'time')
    pr_gpcp_mean = pr_gpcp_season.mean(dim = 'time')

    pr_bias = np.empty((len(pr_cmip6_mean.model),len(pr_cmip6_mean.lat),len(pr_cmip6_mean.lon)))
    pr_bias[:,:,:] = np.nan
    pr_bias = xr.DataArray(pr_bias, coords=pr_cmip6_mean.coords)

    pval = np.empty((len(pr_cmip6_mean.lat),len(pr_cmip6_mean.lon)))
    pval[:,:] = np.nan
    ds_coords = pr_cmip6_mean[0,:,:]
    pval = xr.DataArray(pval, coords=ds_coords.coords)
    pval = pval.drop('model_name')
    
    pr_bias = pr_cmip6_mean - pr_gpcp_mean

    for i in range(len(pval.lat)):
        for j in range(len(pval.lon)):
            pval[i,j] = np.float64(stats.ttest_1samp(pr_bias[:,i,j],popmean=0).pvalue)

    pr_bias_toplot = pr_bias.mean(dim = 'model')

    # pr_cmip6_mmmean = pr_cmip6_season.mean(dim='model').mean(dim = 'time')
    # pr_gpcp_mean = pr_gpcp_season.mean(dim = 'time')

    # pr_bias = pr_cmip6_mmmean - pr_gpcp_mean

    plt.rcParams['hatch.linewidth'] = 0.05 # To control size of hatching!!!!
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    fig.suptitle(title)
    ax.set_xlim([-70,180])
    ax.set_ylim([-20,40])
    ax.set_title(season_title[season],loc='right',fontsize=9)
    ax.set_title(title_ssp[ssp],loc='left',fontsize=9)
    z1_plot = ax.pcolormesh(pr_bias_toplot.lon,pr_bias_toplot.lat,pr_bias_toplot.values,transform=ccrs.PlateCarree(central_longitude=0), cmap= cm.BrBG) #matplotlib.cm.get_cmap('cmo.tarn'))
    zconf_plot = plt.contourf(pval.lon, pval.lat, pval, levels=[0, 0.05],colors='none', hatches=['...'], transform=ccrs.PlateCarree(central_longitude= 0))

    z1_plot.set_clim((-3,3))
    axins = inset_axes(ax,width="70%",  height="10%",loc='lower center',borderpad=-4)
    cbar = fig.colorbar(z1_plot, cax=axins, shrink=0.8, orientation='horizontal')
    cbar.set_label('CMIP6 precipitation bias [mm/day]',fontsize = 11)
    cbar.ax.tick_params(labelsize = 9)
    ax.coastlines(resolution = '50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude= 0), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9, 'color': 'black'}
    gl.ylabel_style = {'size': 9, 'color': 'black'}

    # # BES
    ax.plot([90,90],[5,20],color = 'red',linewidth='1')
    ax.plot([130,130],[5,20],color = 'red',linewidth='1')
    ax.plot([90,130],[5,5],color = 'red',linewidth='1')
    ax.plot([90,130],[20,20],color = 'red',linewidth='1')

    # BES
    # ax.plot([90,90],[5,20],color = 'black',linewidth='1')
    # ax.plot([130,130],[5,20],color = 'black',linewidth='1')
    # ax.plot([90,130],[5,5],color = 'black',linewidth='1')
    # ax.plot([90,130],[20,20],color = 'black',linewidth='1')

    #fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/CMIP6_precipitation_bias_'+season+'_season.png')
    #fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/pwpt/CMIP6_precipitation_bias_'+season+'_season.png')
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_precipitation_bias_'+season+'_season.png',dpi = 300)



breakpoint()

pr_cmip6_BES = pr_cmip6.where((pr_cmip6.lat > 5) & (pr_cmip6.lat < 20) & (pr_cmip6.lon > 270) & (pr_cmip6.lon < 310), drop=True)
pr_gpcp_BES = pr_gpcp.where((pr_gpcp.lat > 5) & (pr_gpcp.lat < 20) & (pr_gpcp.lon > 270) & (pr_gpcp.lon < 310), drop=True)

pr_cmip6_BES_mean = boxmean(pr_cmip6_BES.mean(dim = 'time'))
pr_gpcp_BES_mean = boxmean(pr_gpcp_BES.mean(dim = 'time'))

models = pr_cmip6_BES.model_name

model_bias = np.abs(pr_cmip6_BES_mean - pr_gpcp_BES_mean)
model_bias_quant = model_bias.quantile((0.25,0.75))

for i in range(0,len(pr_cmip6_BES_mean.model_name)):
    strong_model_bias_idx = np.where(model_bias>=model_bias_quant[1])[0]
    weak_model_bias_idx = np.where(model_bias<=model_bias_quant[0])[0]

print('Models with strong precipitation bias are:')
for i in range(len(strong_model_bias_idx)): 
    print('-', models[int(strong_model_bias_idx[i])].values)

print('Models with weak precipitation bias are:')
for i in range(len(weak_model_bias_idx)): 
    print('-', models[int(weak_model_bias_idx[i])].values)

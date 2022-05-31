# This script computes the multimodel mean for the 1991-2020 period for P-E moist_conv mass_conv moist_adv. These variables have been computed before with scripts in folder: Moisture_analysis

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
#import dask

#------------------------INPUT------------------------#


# SETTINGS PROGRAM

# -------------Things to choose:----------------------
# Variables
var1 = 'moist_conv' # P-E moist_conv mass_conv moist_adv

domain = 'tropics'

sign = True

quiver = True # To add wind vectors to the maps
var_quiver = 'wind850'

contour = False # To include contour in the first loop: 30yrs averages, not change in time
var_contour1 = 'psl'

ssp = 'ssp585'

# ----------------------------------------------------

# list things to do
year_period_list = ['2005','2085'] 
season_list = ['dry','wet']#,'dry','wet']
#which_models_list = ['all']#,'dry','wet','dry-wet']

if domain == 'tropics':
    lat_matplot = [-10,40] 
    lon_matplot = [70,160] 
    vector_regrid = 35
elif domain == 'NAtl_NPac':
    lat_matplot = [-20,60]#[-20,60] 
    lon_matplot = [-70,180] 
    vector_regrid = 15  
elif domain == 'tropical_NAtl_NPac':
    lat_matplot = [-10,30]#[-20,60] 
    lon_matplot = [-80,170] 
    vector_regrid = 15 

latit = {'BES': [5,20], 'NAtl_NPac': [-20,60],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40],'tropical_NAtl_NPac':[-10,30] }
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[120,170],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340],'tropical_NAtl_NPac':[130,350]}

latitude = latit[domain]
longitude = longit[domain]



# Paths variables
paths = {'pr': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/pr/',
        'hfls': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/hfls/',
        'moist_conv': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_conv/',
        'moist_mass_conv':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_mass_conv/',
        'moist_adv': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_adv/'}


# Defining new colors 
top = cm.get_cmap('Oranges_r', 128); bottom = cm.get_cmap('Greens', 128); newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128)))); orange_green = ListedColormap(newcolors, name='OrangeGreen')

# SETTINGS FIGURES

# Titles
title = {'2005': 'CMIP6 mmmean (1991-2020). ','2085': 'CMIP6 mmmean (2071-2100). '}
title_var1 = {'wind_shear': 'Zonal wind shear (200hPa - 850 hPa). ', 'P-E': 'Surface P-E. ','hus850': 'Specific humidity. ','pr': 'Precipitation. ','moist_conv': 'Moisture convergence by the mean flow.  ','tas': 'Temperature at 2m. ','wap500': 'Omega at 500 hPa. ','mass_conv': 'Moisture mass convergence.  ', 'moist_adv': 'Moisture advection. '}
season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}

# colors
color_diff_year = {'wind_shear':matplotlib.cm.get_cmap('cmo.diff'), 'P-E':orange_green, 'pr':cm.BrBG,'tas':matplotlib.cm.get_cmap('cmo.balance'),'wap500':cm.coolwarm,'moist_conv': orange_green, 'mass_conv': orange_green, 'moist_adv': orange_green}
# limits
c_lim_diff_year = {'P-E':(-5,5),'moist_conv':(-5,5),'mass_conv':(-5,5),'moist_adv': (-5,5)} 
c_label = {'P-E': r'Scaled $\Delta$ Precipitation - Evaporation [mm/day $^{\circ}C$]', 'moist_conv': r'Scaled $\Delta$ Moisture convergence [mm/day $^{\circ}C$]', 'mass_conv': r'Scaled $\Delta$ Mass convergence [mm/day $^{\circ}C$]', 'moist_adv': r'Scaled $\Delta$ Moisture advection [mm/day $^{\circ}C$]'}
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
  
    if 'plev_bnds' in ds:
        ds =  ds.drop('plev_bnds')
    if 'plev' in ds:
        ds = ds.drop('plev')
    if 'height' in ds:
        ds = ds.drop('height')
    if 'file_qf' in ds:
        ds = ds.drop('file_qf')
    if 'lon_bnds' in ds:
        ds = ds.drop('lon_bnds')
    if 'lat_bnds' in ds: 
        ds = ds.drop('lat_bnds')
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


for season in ['dry','wet']:

    # Import variables

    # Variable 1: Colormap variable
    if var1 == 'P-E':
        pr_2005 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2085 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2005 = pr_2005.load(); pr_2085 = pr_2085.load()

        hfls_2005 = xr.open_mfdataset(paths['hfls']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        hfls_2085 = xr.open_mfdataset(paths['hfls']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        hfls_2005 = hfls_2005.load(); hfls_2085 = hfls_2085.load()

        ds_2005 = pr_2005.pr*86400 - hfls_2005.hfls * 86400/(1000*2455) # P-E in mm/day 
        ds_2085 = pr_2085.pr*86400 - hfls_2085.hfls * 86400/(1000*2455)

    elif var1 == 'moist_conv':
        moist_conv_2005 = xr.open_mfdataset(paths['moist_conv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        moist_conv_2085 = xr.open_mfdataset(paths['moist_conv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        moist_conv_2005 = moist_conv_2005.load(); moist_conv_2085 = moist_conv_2085.load()
        ds_2005 = moist_conv_2005.conv[:,0,:,:]
        ds_2085 = moist_conv_2085.conv[:,0,:,:]

    elif var1 == 'mass_conv':
        mass_conv_2005 = xr.open_mfdataset(paths['moist_mass_conv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        mass_conv_2085 = xr.open_mfdataset(paths['moist_mass_conv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        mass_conv_2005 = mass_conv_2005.load(); mass_conv_2085 = mass_conv_2085.load()
        ds_2005 = mass_conv_2005.mass_conv[:,0,:,:]
        ds_2085 = mass_conv_2085.mass_conv[:,0,:,:]

    elif var1 == 'moist_adv':
        moist_adv_2005 = xr.open_mfdataset(paths['moist_adv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        moist_adv_2085 = xr.open_mfdataset(paths['moist_adv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        moist_adv_2005 = moist_adv_2005.load(); moist_adv_2085 = moist_adv_2085.load()
        ds_2005 = moist_adv_2005.adv[:,0,:,:]
        ds_2085 = moist_adv_2085.adv[:,0,:,:]



    
    ds_2005 = ds_2005.sortby(ds_2005.model_name)

    # Var 1
    pval = np.empty((len(ds_2005.lat),len(ds_2005.lon))) #ds_2005[0,:,:]
    pval[:,:] = np.nan
    ds_coord = ds_2005[0,:,:]
    pval = xr.DataArray(pval, coords=ds_coord.coords)
    pval = pval.drop('model_name')
    #breakpoint()
    for i in range(len(pval.lat)):
        for j in range(len(pval.lon)):
            pval[i,j] = np.float64(stats.ttest_1samp(ds_2005[:,i,j],popmean=0).pvalue)

    #breakpoint()

    ds_toplot = ds_2005.mean(dim = 'model', skipna = False)

       
    #Figure
    plt.rcParams['hatch.linewidth'] = 0.05 # To control size of hatching!!!!
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    fig.suptitle('CMIP6 mmm 2005 (1991-2020). '+title_var1[var1]) 
    ax.set_title(season_title[season],loc='right',fontsize=14)
    ax.set_title(title_ssp[ssp],loc='left',fontsize=14)

    z1_plot = ax.pcolormesh(ds_toplot.lon, ds_toplot.lat, ds_toplot,transform=ccrs.PlateCarree(central_longitude= 0), cmap= color_diff_year[var1])
    if sign == True: 
        zconf_plot = plt.contourf(pval.lon, pval.lat, pval, levels=[0, 0.05],colors='none', hatches=['....'], transform=ccrs.PlateCarree(central_longitude= 0))

    ax.coastlines(resolution = '50m')

    # BES
    ax.plot([90,90],[5,20],color = 'red',linewidth='1')
    ax.plot([130,130],[5,20],color = 'red',linewidth='1')
    ax.plot([90,130],[5,5],color = 'red',linewidth='1')
    ax.plot([90,130],[20,20],color = 'red',linewidth='1')

    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude= 0), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 11  , 'color': 'black'}
    gl.ylabel_style = {'size': 11, 'color': 'black'}

    #ax.gridlines(draw_labels=True)
    ax.set_xlim(lon_matplot)
    ax.set_ylim(lat_matplot)
    z1_plot.set_clim(c_lim_diff_year[var1])
    axins = inset_axes(ax,width="70%",  height="8%",loc='lower center',borderpad=-4.5)
    
    cbar = fig.colorbar(z1_plot, cax=axins, shrink=0.5, orientation='horizontal')
    cbar.set_label(c_label[var1],fontsize=14)
    cbar.ax.tick_params(rotation=45,labelsize = 11)

    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_spatialplot_'+var1+'_2005_'+ season+'_season.png',dpi = 1200)
    print('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_spatialplot_'+var1+'_2005_'+ season+'_season.png')
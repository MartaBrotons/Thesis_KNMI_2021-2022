# This script computes the scaled change in tas, hus, hur and wap averaged over the Caribbean region for the whole atmospheric column (Figures 18 and 19)

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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

#------------------------INPUT------------------------#

# SETTINGS PROGRAM

# -------------Things to choose:----------------------
# Variables

var = 'hur' # hur hus ta wap
n = 58986
ssp = 'ssp585'

# ----------------------------------------------------

# list things to do
season_list = ['annual','dry','wet']
number_list = [0,1,2]
year_period_list = ['2005', '2085']

# Paths variables
paths = {'hus':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/hus/*'+ssp+'*.nc' ,
        'ta':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/ta/*'+ssp+'*.nc',
        'wap': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/wap/*'+ssp+'*.nc',
        'hur': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/hur/*'+ssp+'*.nc'}

path_area_mean = {'tas_40N_40S': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*'+ssp+'*.nc'}

# SETTINGS FIGURES

period_date = {'2005': ['1-1-1991','31-12-2020'], '2085' : ['1-1-2071','31-12-2100']}


xlabel = {'hus': '$\Delta$ Specific humidity', 'ta': r'$\Delta$ Temperature [$^{\circ}$C/$^{\circ}$C]', 'wap': '$\Delta$ Omega [Pa/s $^{\circ}$C]', 'hur': '$\Delta$ Relative humidity [%/ $^{\circ}$C]'}

# Titles

title_var = {'hus': 'Specific humidity','ta': 'Air temperature', 'wap': 'Vertical motion','hur': 'Relative humidity'}
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

def preprocessing_vertical(ds):
    """
    Function to call when importing netcdf files with xr.open_mfdataset.
    - Sets all the datasets to the same datetime and time axis
    - Deletes from the files inconsistent dimensions/coordinates and saves models name
    - Saves the name of the model and the rpf number
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

tas_40 = xr.open_mfdataset(path_area_mean['tas_40N_40S'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
tas_40 = tas_40.load()
tas_40 = tas_40.tas
tas_40 = tas_40.sortby(tas_40.model_name)

ds = xr.open_mfdataset(paths[var],combine = 'nested', concat_dim = 'model',preprocess = preprocessing_vertical)
ds = ds.load()
ds = ds.sortby(ds.model_name)

if var == 'ta': ds = ds.ta[:,:,:,0,0]-273.15; color_line =  'darkred'; color_range = 'lightcoral'; x_0 = 1; x_limit = (-3.5,2.5)
elif var == 'hus': ds = ds.hus[:,:,:,0,0]; color_line = 'steelblue'; color_range = 'lightskyblue'; x_0 = 0
elif var == 'wap': ds = ds.wap[:,:,:,0,0]; color_line = 'chocolate'; color_range =   'sandybrown'; x_0 = 0
elif var == 'hur': ds = ds.hur[:,:,:,0,0]; color_line ='steelblue'; color_range = 'lightskyblue'; x_0 = 0

plevs = [100000.,  92500.,  85000.,  70000.,  60000.,  50000.,  40000.,  30000., 25000.,  20000.,  15000.,  10000.,   7000.,   5000.,   3000.,   2000., 1000.,    500.,    100.]
ds['plev'] = plevs

#----------------------------------#

# Start computations

breakpoint()

for season,number in zip(season_list,number_list):

    ds_season = time_period(ds,season)

    #fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    plt.figure(n+number)
    


    ds_2005 = ds_season.sel(time=slice(period_date['2005'][0], period_date['2005'][1])).mean(dim = 'time') 
    ds_2085 = ds_season.sel(time=slice(period_date['2085'][0], period_date['2085'][1])).mean(dim = 'time') 
    if (var == 'hus') or (var == 'hur'): ds_2005 = ds_2005.where((ds_2005 >=0), drop = True); ds_2085 = ds_2085.where((ds_2085 >=0), drop = True)
    tas_40_2005 = tas_40.sel(time=slice(period_date['2005'][0], period_date['2005'][1])).mean(dim = 'time') 
    tas_40_2085 = tas_40.sel(time=slice(period_date['2085'][0], period_date['2085'][1])).mean(dim = 'time') 
    ds_toplot = (ds_2085.mean(dim = 'model') - ds_2005.mean(dim = 'model'))/(tas_40_2085.mean(dim = 'model') - tas_40_2005.mean(dim = 'model'))
    ds_toquant = (ds_2085 - ds_2005)/(tas_40_2085 - tas_40_2005)
    

    
    plt.gca().invert_yaxis()
    plt.plot(ds_toplot,ds_toplot.plev*1e-2,color = color_line )
    plt.fill_betweenx(ds_toplot.plev*1e-2,ds_toquant.quantile(0.05,dim='model'),ds_toquant.quantile(0.95,dim='model'),color=color_range,linewidth=0)
    plt.axvline(x= x_0, color='black', ls = '--', lw = 0.5)
    # plt.axvline(x= 0, color='grey', ls = '--', lw = 0.5)
    # plt.xlim(x_limit)
    plt.gca().invert_yaxis()

    plt.xlabel(xlabel[var],fontsize=12)
    plt.ylabel('Pressure level [hPa]',fontsize=12)
    plt.title(season_title[season],loc='right',fontsize=11)
    plt.title(title_ssp[ssp],loc='left',fontsize=11)

    # fig.suptitle('30 years mean vertical profile. '+title_var[var])


    percent =  mpatches.Patch( color=color_range, label=r'$90\%$ CMIP6')
    line = mlines.Line2D([], [], color=color_line, label=r'CMIP6 mmmean')
    plt.legend(handles=[percent,line])
                    
    plt.gca().invert_yaxis()
    # plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/CMIP6_vertical_profile_2085-2005_'+var+'_'+season+'_season_'+ssp+'.png')
    # print('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/CMIP6_vertical_profile_2085-2005_'+var+'_'+season+'_season_'+ssp+'.png')

    plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_vertical_profile_2085-2005_'+var+'_'+season+'_season_'+ssp+'.png',dpi = 300)
    print('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_vertical_profile_2085-2005_'+var+'_'+season+'_season_'+ssp+'.png')
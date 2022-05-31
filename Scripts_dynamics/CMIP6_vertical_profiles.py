# This script computes 1991-2020 and 2071-2100 averaged for tas, hus, hur and wap averaged over the Caribbean region for the whole atmospheric column (Figure A5)
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

var = 'tas' # hur hus ta wap
n = 4567890098765678
ssp = 'ssp585'

# ----------------------------------------------------

# list things to do
season_list = ['annual']#,'dry','wet']
number_list = [0,1,2]
year_period_list = ['2005', '2085']

# Paths variables
paths = {'hus':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/hus/*'+ssp+'*.nc' ,
        'ta':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/ta/*'+ssp+'*.nc',
        'wap': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/wap/*'+ssp+'*.nc',
        'hur': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6_vertical/00_join_CMIP6_vertical_data/hur/*'+ssp+'*.nc'}

# SETTINGS FIGURES

period_date = {'2005': ['1-1-1991','31-12-2020'], '2085' : ['1-1-2071','31-12-2100']}


xlabel = {'hus': 'Specific humidity', 'ta': r'Temperature [$^{\circ}$C]', 'wap': 'Omega [Pa/s]', 'hur': 'Relative humidity [%]'}
#xlim = {'hus':(,),'tas':(,),'wap':''}

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

#------------------------CALCULATIONS------------------------#

# Import variables
breakpoint()

ds = xr.open_mfdataset(paths[var],combine = 'nested', concat_dim = 'model',preprocess = preprocessing_vertical)
ds = ds.load()

if var == 'ta': ds = ds.ta[:,:,:,0,0]-273.15; color_line = {'2005': 'black', '2085': 'darkred'}; color_range = {'2005': 'darkgrey', '2085': 'lightcoral'}
elif var == 'hus': ds = ds.hus[:,:,:,0,0]; color_line = {'2005': 'black', '2085': 'steelblue'}; color_range = {'2005': 'darkgrey', '2085': 'lightskyblue'}
elif var == 'wap': ds = ds.wap[:,:,:,0,0]; color_line = {'2005': 'black', '2085': 'peru'}; color_range = {'2005': 'darkgrey', '2085': 'sandybrown'}
elif var == 'hur': ds = ds.hur[:,:,:,0,0]; color_line = {'2005': 'black', '2085': 'steelblue'}; color_range = {'2005': 'darkgrey', '2085': 'lightskyblue'}

plevs = [100000.,  92500.,  85000.,  70000.,  60000.,  50000.,  40000.,  30000., 25000.,  20000.,  15000.,  10000.,   7000.,   5000.,   3000.,   2000., 1000.,    500.,    100.]
ds['plev'] = plevs

#----------------------------------#

# Start computations

for season,number in zip(season_list,number_list):

    ds_season = time_period(ds,season)

    #fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    plt.figure(n+number)
    

    for year_period in year_period_list: 

        ds_period = ds_season.sel(time=slice(period_date[year_period][0], period_date[year_period][1])).mean(dim = 'time') 
        if var == 'hus': ds_period = ds_period.where((ds_period >=0), drop = True)
        ds_toplot = ds_period.mean(dim = 'model')

        plt.gca().invert_yaxis()
        plt.plot(ds_toplot,ds_toplot.plev*1e-2,color = color_line[year_period] )
        plt.fill_betweenx(ds_toplot.plev*1e-2,ds_period.quantile(0.05,dim='model'),ds_period.quantile(0.95,dim='model'),color=color_range[year_period],linewidth=0)
        plt.gca().invert_yaxis()

    plt.xlabel(xlabel[var])
    plt.ylabel('Pressure level [hPa]')
    plt.title(season_title[season],loc='right',fontsize=9)
    plt.title(title_ssp[ssp],loc='left',fontsize=9)
    # fig.suptitle('30 years mean vertical profile. '+title_var[var])

    percent_2005 =  mpatches.Patch( color=color_range['2005'], label=r'$90\%$ CMIP6 for 2005 ')
    percent_2085 =  mpatches.Patch( color=color_range['2085'], label=r'$90\%$ CMIP6 for 2085 ')
    line_2005 = mlines.Line2D([], [], color=color_line['2005'], label=r'CMIP6 mmmean for 2005')
    line_2085 = mlines.Line2D([], [], color=color_line['2085'], label=r'CMIP6 mmmean for 2085')
    plt.legend(handles=[percent_2005,percent_2085,line_2005,line_2085])
                
    plt.gca().invert_yaxis()
    plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_vertical_profile_'+var+'_'+season+'_season_'+ssp+'.png',dpi = 300)
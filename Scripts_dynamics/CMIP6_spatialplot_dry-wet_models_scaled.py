# This script computes the dry models multimodel mean - wet models multimodel mean. Dry and wet models list come from CMIP6_dry_wet_models_definition.py
# Variables: tas pr wap500 P-E moist_conv mass_conv moist_adv wind_shear


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

# -------------Things to choose:----------------------
# Variables
var = 'pr' # tas pr wap500 P-E moist_conv mass_conv moist_adv wind_shear
ssp = 'ssp585'
domain = 'NAtl_NPac'
sign = True

quiver = False # To add wind vectors to the maps
var_quiver = 'wind850'
# ----------------------------------------------------

if domain == 'tropics':
    lat_matplot = [-10,40] 
    lon_matplot = [70,160] 
    vector_regrid = 25
elif domain == 'NAtl_NPac':
    lat_matplot = [-20,60]#[-20,60] 
    lon_matplot = [-70,180] 
    vector_regrid = 30  

latit = {'BES': [5,20], 'NAtl_NPac': [-20,60],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40]}
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[120,170],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340]}

# SETTINGS PROGRAM

name_quiv = ''
title_quiv = ''
length_vector = 1
label_vector = r'1 m/s $^{\circ}$C'
cent_lon = 0 #always 0 I think, if not also 180
loc_quiverkey = [0.04,-0.25] 


# 30 and 70 percentiles

wet_models = ['CIESM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3','FGOALS-f3-L','INM-CM4-8','INM-CM5-0','NESM3']
dry_models = ['ACCESS-ESM1-5','AWI-CM-1-1-MR','CESM2','CanESM5','CanESM5-CanOE','FGOALS-g3','IPSL-CM6A-LR','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0']

latitude = latit[domain]
longitude = longit[domain]

dry_wet_models = {'wet':wet_models,'dry':dry_models}

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
        'moist_adv': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/moist_adv/',
        'psl': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/psl/' }

path_area_mean = {'tas_40N_40S': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*'+ssp+'*.nc'}


# Defining new colors 
top = cm.get_cmap('Blues_r', 128); bottom = cm.get_cmap('Oranges', 128);newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128)))); blue_orange = ListedColormap(newcolors, name='BlueOrange')
top = cm.get_cmap('Oranges_r', 128); bottom = cm.get_cmap('Greens', 128); newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128)))); orange_green = ListedColormap(newcolors, name='OrangeGreen')
top = cm.get_cmap('YlOrBr_r', 128); bottom = cm.get_cmap('Blues', 128); newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128)))); yellow_blue = ListedColormap(newcolors, name='YellowBlue')

# SETTINGS FIGURES

# Titles
season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
title_models = {'dry': 'Dry models.  Difference ref. period (1991-2020) and 2085 (2071-2100)', 'wet':'Wet models.  Difference ref. period (1991-2020) and 2085 (2071-2100)'  }
title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}

title_var = {'wind_shear': 'Zonal wind shear (200hPa - 850 hPa). ', 'P-E': 'Surface P-E. ','hus850': 'Specific humidity. ','pr': 'Precipitation. ','moist_conv': 'Moisture convergence by the mean flow.  ','tas': 'Temperature at 2m. ','wap500': 'Omega at 500 hPa. ','mass_conv': 'Moisture mass convergence.  ', 'moist_adv': 'Moisture advection. ', 'psl': 'Mean sea level pressure '}
title_quiver = {'wind850': r'Scaled $\Delta$ wind at 850 hPa. ', 'wind200': r'Scaled $\Delta$ wind at 200 hPa. ' }
# colors
color_diff_year = {'wind_shear':matplotlib.cm.get_cmap('cmo.diff'), 'P-E':yellow_blue, 'pr':cm.BrBG,'tas':matplotlib.cm.get_cmap('cmo.balance'),'wap500':cm.coolwarm,'moist_conv': yellow_blue, 'mass_conv': yellow_blue, 'moist_adv': yellow_blue, 'psl': blue_orange}
# limits
c_lim_diff_year = {'wind_shear':(-0.7,0.7),'P-E':(-.4,.4),'pr': (-0.7,0.7),'tas':(-0.5,0.5),'wap500':(-0.005, 0.005),'moist_conv':(-.4,.4),'mass_conv':(-.4,.4),'moist_adv': (-.4,.4),'psl':(-0.6,0.6)} #moist_conv850'(-0.05,0.05)} for 225dP # colorbar limits 2085-2005t models plot #'850': (-2,2), '200': (-4,4),
c_label = {'wind_shear': r'Scaled $\Delta$ VWS [m/s $^{\circ}C$]', 'P-E': r'Scaled $\Delta$ Precipitation - Evaporation [mm/day $^{\circ}C$]', 'pr': r'Scaled $\Delta$ Precipitation [mm/day $^{\circ}C$]','wap500': r'Scaled $\Delta$ Omega [Pa/s $^{\circ}C$]','tas': r'Scaled $\Delta$ Temperature at 2m [$^{\circ}C$/$^{\circ}C$]', 'moist_conv': r'Scaled $\Delta$ Moisture convergence [mm/day $^{\circ}C$]', 'mass_conv': r'Scaled $\Delta$ Mass convergence [mm/day $^{\circ}C$]', 'moist_adv': r'Scaled $\Delta$ Moisture advection [mm/day $^{\circ}C$]', 'psl': r'Scaled $\Delta$ Mean sea level pressure [hPa/ $^{\circ}C$]'}

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
    return ds

def preprocessing_lat_lon_moist(ds):
    """
    Function to call when importing netcdf files with xr.open_mfdataset.
    - Sets all the datasets to the same datetime and time axis
    - Selects the domain we are interested on
    - Deletes from the files inconsistent dimensions/coordinates and saves models name
    """
    if ds.indexes['time'].inferred_type != "datetime64":
            ds['time'] = ds.indexes['time'].to_datetimeindex()
    ds['time'] = pd.date_range(start='1950-01-01',end='2100-12-16',periods=len(ds.time))
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

# Import variable


tas_40 = xr.open_mfdataset(path_area_mean['tas_40N_40S'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
tas_40 = tas_40.load()
tas_40 = tas_40.tas
tas_40 = tas_40.sortby(tas_40.model_name)


for season in ['dry','wet']:


    # Variable 1: Colormap variable
    if var == 'wind_shear':
        ua850_2005 = xr.open_mfdataset(paths['ua850']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        ua850_2085 = xr.open_mfdataset(paths['ua850']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        ua850_2005 = ua850_2005.load(); ua850_2085 = ua850_2085.load()

        ua200_2005 = xr.open_mfdataset(paths['ua200']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        ua200_2085 = xr.open_mfdataset(paths['ua200']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        ua200_2005 = ua200_2005.load(); ua200_2085 = ua200_2085.load()

        ds_2005 = ua200_2005.ua[:,0,:,:] - ua850_2005.ua[:,0,:,:]
        ds_2085 = ua200_2085.ua[:,0,:,:] - ua850_2085.ua[:,0,:,:]

    elif var == 'P-E':
        pr_2005 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2085 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2005 = pr_2005.load(); pr_2085 = pr_2085.load()

        hfls_2005 = xr.open_mfdataset(paths['hfls']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        hfls_2085 = xr.open_mfdataset(paths['hfls']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        hfls_2005 = hfls_2005.load(); hfls_2085 = hfls_2085.load()

        ds_2005 = pr_2005.pr*86400 - hfls_2005.hfls * 86400/(1000*2455) # P-E in mm/day 
        ds_2085 = pr_2085.pr*86400 - hfls_2085.hfls * 86400/(1000*2455)

    elif var == 'pr':
        pr_2005 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2085 = xr.open_mfdataset(paths['pr']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        pr_2005 = pr_2005.load(); pr_2085 = pr_2085.load()
        ds_2005 = pr_2005.pr*86400
        ds_2085 = pr_2085.pr*86400


    elif var == 'tas':
        tas_2005 = xr.open_mfdataset(paths['tas']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        tas_2085 = xr.open_mfdataset(paths['tas']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        tas_2005 = tas_2005.load(); tas_2085 = tas_2085.load()
        ds_2005 = tas_2005.tas
        ds_2085 = tas_2085.tas

    elif var == 'wap500':
        wap500_2005 = xr.open_mfdataset(paths['wap500']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        wap500_2085 = xr.open_mfdataset(paths['wap500']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        wap500_2005 = wap500_2005.load(); wap500_2085 = wap500_2085.load()
        ds_2005 = wap500_2005.wap[:,0,:,:]
        ds_2085 = wap500_2085.wap[:,0,:,:]


    elif var == 'moist_conv':
        moist_conv_2005 = xr.open_mfdataset(paths['moist_conv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        moist_conv_2085 = xr.open_mfdataset(paths['moist_conv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        moist_conv_2005 = moist_conv_2005.load(); moist_conv_2085 = moist_conv_2085.load()
        ds_2005 = moist_conv_2005.conv[:,0,:,:]
        ds_2085 = moist_conv_2085.conv[:,0,:,:]

    elif var == 'mass_conv':
        mass_conv_2005 = xr.open_mfdataset(paths['moist_mass_conv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        mass_conv_2085 = xr.open_mfdataset(paths['moist_mass_conv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        mass_conv_2005 = mass_conv_2005.load(); mass_conv_2085 = mass_conv_2085.load()
        ds_2005 = mass_conv_2005.mass_conv[:,0,:,:]
        ds_2085 = mass_conv_2085.mass_conv[:,0,:,:]

    elif var == 'moist_adv':
        moist_adv_2005 = xr.open_mfdataset(paths['moist_adv']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
        moist_adv_2085 = xr.open_mfdataset(paths['moist_adv']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
        moist_adv_2005 = moist_adv_2005.load(); moist_adv_2085 = moist_adv_2085.load()
        ds_2005 = moist_adv_2005.adv[:,0,:,:]
        ds_2085 = moist_adv_2085.adv[:,0,:,:]

    elif var == 'psl':
        psl_2005 = xr.open_mfdataset(paths['psl']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        psl_2085 = xr.open_mfdataset(paths['psl']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon) 
        psl_2005 = psl_2005.load(); psl_2085 = psl_2085.load()
        ds_2005 = psl_2005.psl*1e-2
        ds_2085 = psl_2085.psl*1e-2

    # tas to scale
    ds_2005 = ds_2005.sortby(ds_2005.model_name)
    ds_2085 = ds_2085.sortby(ds_2085.model_name)
    tas_40_2005 = tas_40.where((tas_40.time.dt.year>=1991) & (tas_40.time.dt.year<=2020),drop=True).mean(dim = 'time')
    tas_40_2085 = tas_40.where((tas_40.time.dt.year>=2071) & (tas_40.time.dt.year<=2100),drop=True).mean(dim = 'time')

    # Variable 2: Quiver variable
    if quiver == True: 
        title_quiv = title_quiver[var_quiver]
        name_quiv = '_'+var_quiver
        if var_quiver == 'wind850':
            ua850_2005 = xr.open_mfdataset(paths['ua850']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
            ua850_2085 = xr.open_mfdataset(paths['ua850']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
            ua850_2005 = ua850_2005.load(); ua850_2085 = ua850_2085.load()

            va850_2005 = xr.open_mfdataset(paths['va850']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
            va850_2085 = xr.open_mfdataset(paths['va850']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
            va850_2005 = va850_2005.load(); va850_2085 = va850_2085.load()

            ua_2005 = ua850_2005.ua[:,0,:,:]; ua_2085 =  ua850_2085.ua[:,0,:,:]
            va_2005 = va850_2005.va[:,0,:,:]; va_2085 =  va850_2085.va[:,0,:,:]

        if var_quiver == 'wind200':
            ua200_2005 = xr.open_mfdataset(paths['ua200']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
            ua200_2085 = xr.open_mfdataset(paths['ua200']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
            ua200_2005 = ua200_2005.load(); ua200_2085 = ua200_2085.load()

            va200_2005 = xr.open_mfdataset(paths['va200']+'*_'+season+'_season_2005.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)#preprocessing_save_model_names)
            va200_2085 = xr.open_mfdataset(paths['va200']+'*_'+season+'_season_2085.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing_lat_lon)
            va200_2005 = va200_2005.load(); va200_2085 = va200_2085.load()

            ua_2005 = ua200_2005.ua[:,0,:,:]; ua_2085 =  ua200_2085.ua[:,0,:,:]
            va_2005 = va200_2005.va[:,0,:,:]; va_2085 =  va200_2085.va[:,0,:,:]

        ua = (ua_2085.sortby(ua_2085.model_name)- ua_2005.sortby(ua_2005.model_name))/(tas_40_2085-tas_40_2005)
        va = (va_2085.sortby(va_2085.model_name)- va_2005.sortby(va_2005.model_name))/(tas_40_2085-tas_40_2005)


    # Var 1
    ds = (ds_2085 - ds_2005)/(tas_40_2085-tas_40_2005)

    # Select dry and wet models

    models_wet_idx = np.zeros(len(wet_models))
    for i in range(0,len(wet_models)):
        models_wet_idx[i] = np.where(ds.model_name == wet_models[i])[0]

    models_dry_idx = np.zeros(len(dry_models))
    for i in range(0,len(dry_models)):
        models_dry_idx[i] = np.where(ds.model_name == dry_models[i])[0]

    models_wet_idx = models_wet_idx.astype(int)
    models_dry_idx = models_dry_idx.astype(int)

    ds_dry_models = ds.sel(model = models_dry_idx)
    ds_wet_models = ds.sel(model = models_wet_idx)

    if quiver == True:
        ua_dry_models = ua.sel(model = models_dry_idx); ua_wet_models = ua.sel(model = models_wet_idx)
        va_dry_models = va.sel(model = models_dry_idx); va_wet_models = va.sel(model = models_wet_idx)
        ua_toplot = ua_dry_models.mean(dim = 'model') - ua_wet_models.mean(dim = 'model')
        va_toplot = va_dry_models.mean(dim = 'model') - va_wet_models.mean(dim = 'model')


    #breakpoint()
    pval_dry = np.empty((len(ds_dry_models.lat),len(ds_dry_models.lon))) #ds_toplot[0,:,:]
    pval_dry[:,:] = np.nan
    ds_coords = ds_dry_models[0,:,:]
    pval_dry = xr.DataArray(pval_dry, coords=ds_coords.coords)
    pval_dry = pval_dry.drop('model_name')

    pval_wet = np.empty((len(ds_wet_models.lat),len(ds_wet_models.lon))) #ds_toplot[0,:,:]
    pval_wet[:,:] = np.nan
    pval_wet = xr.DataArray(pval_wet, coords=ds_coords.coords)
    pval_wet = pval_wet.drop('model_name')

    for i in range(len(pval_dry.lat)):
        for j in range(len(pval_dry.lon)):
            pval_dry[i,j] = np.float64(stats.ttest_1samp(ds_dry_models[:,i,j],popmean=0).pvalue)
            pval_wet[i,j] = np.float64(stats.ttest_1samp(ds_wet_models[:,i,j],popmean=0).pvalue)

    ds_toplot = ds_dry_models.mean(dim = 'model', skipna = False)-ds_wet_models.mean(dim = 'model', skipna = False)
    pval = (pval_dry + pval_wet)/2


    #breakpoint()
    #Figure
    plt.rcParams['hatch.linewidth'] = 0.05 # To control size of hatching!!!!
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    fig.suptitle('Dry - Wet CMIP6. Change 2085 (2071-2100)- 2005 (1991-2020). '+title_var[var]+title_quiv) 
    ax.set_title(season_title[season],loc='right',fontsize=14)
    ax.set_title(title_ssp[ssp],loc='left',fontsize=14)

    z1_plot = ax.pcolormesh(ds_toplot.lon, ds_toplot.lat, ds_toplot,transform=ccrs.PlateCarree(central_longitude= 0), cmap= color_diff_year[var])
    if sign == True: 
        zconf_plot = plt.contourf(pval.lon, pval.lat, pval, levels=[0, 0.05],colors='none', hatches=['....'], transform=ccrs.PlateCarree(central_longitude= 0))

    if quiver == True: 
        z2_plot = ax.quiver(ua_toplot.lon, ua_toplot.lat, ua_toplot.values, va_toplot.values,width=0.002,transform=ccrs.PlateCarree(central_longitude=cent_lon),regrid_shape=vector_regrid)
        ax.quiverkey(z2_plot,loc_quiverkey[0]+0.005,loc_quiverkey[1],length_vector, label_vector)


    ax.coastlines(resolution = '50m')

    # BES
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

    ax.set_xlim(lon_matplot)
    ax.set_ylim(lat_matplot)
    z1_plot.set_clim(c_lim_diff_year[var])
    axins = inset_axes(ax,width="70%",  height="8%",loc='lower center',borderpad=-4.5)
    cbar = fig.colorbar(z1_plot, cax=axins, shrink=0.5, orientation='horizontal')
    cbar.set_label(c_label[var],fontsize=14)
    cbar.ax.tick_params(rotation=45,labelsize = 11)
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_spatialplot_dry-wet_'+var+name_quiv+'_2085-2005_'+ season+'_season.png',dpi = 300)
    print('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/CMIP6_spatialplot_dry-wet_'+var+name_quiv+'_2085-2005_'+ season+'_season.png')
    # fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/pwpt/CMIP6_spatialplot_dry-wet_'+var+name_quiv+'_2085-2005_'+ season+'_season.png')
    # print('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/pwpt/CMIP6_spatialplot_dry-wet_'+var+name_quiv+'_2085-2005_'+ season+'_season.png')

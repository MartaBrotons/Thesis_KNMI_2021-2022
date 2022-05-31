# Python Script to compute 30 years trends for CMIP6 models for the end of the 21st century. Figures not included in the report
# There are two posibbilities: Niño 3.4 region SST trends, and AMOC proxy (NASPG region tas)


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

# --------------Things to choose --------------------
var = 'Niño34' # AMOC Niño34
ssp = 'ssp585'#,'ssp245','ssp126']
# trend = True
scaling = False
# ---------------------------------------------------

# SETTINGS PROGRAM

# List of things to do

season_list = ['annual','dry','wet']
domain = 'Nino34' # BES NAtl_NPac 40N_40S Nino34 Nino12 NASPG tropics

# Paths
paths = {'ts_Nino34_area_mean':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/Nino34/ts/*'+ssp+'*.nc',
        'tas_40N_40S_area_mean': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/40N_40S/tas/*'+ssp+'*.nc',
        'tas_NH_poles_nonNASPG_area_mean':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/NH_poles_nonNASPG/tas/*'+ssp+'*.nc',
        'tas_NASPG_area_mean':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/NASPG/tas/*'+ssp+'*.nc'}

paths_era5 = {'sst':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_SST_glob.nc', #'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_SST_nino34.nc',#
        'tas':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_t2m_glob.nc'}

# definition coordinates
latit = {'BES': [5,20], 'NAtl_NPac': [-20,60],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40]}
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[190,240],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340]}

latitude = latit[domain]
longitude = longit[domain]

marker = ['v','*','D','>','1','s','p','h']

title_ssp = {'ssp126':'SSP1-2.6', 'ssp245': 'SSP2-4.5', 'ssp585': 'SSP5-8.5'}
season_title = {'annual':'Annual','dry': 'Dry season (December-April)','wet':'Wet season (May-November)'}
clabel_trend = {'Niño34': 'Niño 3.4 SST Trend [K/year]', 'AMOC': 'AMOC tas Trend [K/year]' }
clabel_trend_scaled = {'Niño34': 'Scaled Niño 3.4 SST Trend [K/K year]', 'AMOC': 'Scaled AMOC tas Trend [K/K year]' }


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

def boxmean_era(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.latitude))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('latitude','longitude'))#.dropna(dim='latitude').dropna(dim='longitude')
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
breakpoint()
# Input variables

if var == 'Niño34':
    # CMIP6 data:
    ts = xr.open_mfdataset(paths['ts_Nino34_area_mean'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
    ts = ts.load()
    ts = ts.ts
    ts_ref = ts.sel(time=slice('01-01-1991', '31-12-2020')).mean(dim = 'time')
    ts = ts.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    # ts = ts.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    # ts_ref = ts.sel(time=slice('01-01-2071', '31-12-2100')).mean(dim = 'time')
    ds = ts.sortby(ts.model_name)-ts_ref.sortby(ts_ref.model_name)


elif var == 'AMOC':
    # CMIP6 data:
    #tas_NAtl = xr.open_mfdataset(paths['tas_NH_poles_nonNASPG_area_mean'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
    tas_NASPG = xr.open_mfdataset(paths['tas_NASPG_area_mean'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
    #tas_NAtl = tas_NAtl.load()
    #tas_NAtl = tas_NAtl.tas
    tas_NASPG = tas_NASPG.load()
    tas_NASPG = tas_NASPG.tas
    #tas_NAtl_reg = tas_NAtl.sel(time=slice('1-1-1979', '31-12-2020'),drop = True)
    tas_NASPG_ref = tas_NASPG.sel(time=slice('1-1-1991', '31-12-2020'),drop = True).mean(dim = 'time')
    tas_NASPG = tas_NASPG.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    # tas_NAtl_reg = tas_NAtl.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    # tas_NASPG_reg = tas_NASPG.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    ds = tas_NASPG.sortby(tas_NASPG.model_name)-tas_NASPG_ref.sortby(tas_NASPG_ref.model_name) #- tas_NAtl_reg.sortby(tas_NAtl_reg.model_name)

breakpoint()
ds = ds.sortby(ds.model_name)

if scaling == True:
    sc = '_sc_pr'
    title_sc = 'Scaled Precipitation.'
    title_var1_sc = 'Scaled '
    # 30 and 70 percentiles
    wet_models = ['CIESM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3','FGOALS-f3-L','INM-CM4-8','INM-CM5-0','NESM3']
    dry_models = ['ACCESS-ESM1-5','AWI-CM-1-1-MR','CESM2','CanESM5','CanESM5-CanOE','FGOALS-g3','IPSL-CM6A-LR','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0']
    clabel_trend = clabel_trend_scaled

    # CMIP6
    tas_40N_40S = xr.open_mfdataset(paths['tas_40N_40S_area_mean'],concat_dim = 'model',combine = 'nested', preprocess = preprocessing_area_mean)
    tas_40N_40S = tas_40N_40S.load()
    tas_40N_40S = tas_40N_40S.tas
    tas_40N_40S_reg = tas_40N_40S.sel(time=slice('1-1-2059', '31-12-2100'),drop = True)
    tas_40N_40S_reg['time'] = ds.time 
    tas_40N_40S_reg = tas_40N_40S_reg.sortby(tas_40N_40S_reg.model_name)
    ds = ds/tas_40N_40S_reg



elif scaling == False:
    sc = ''
    title_sc = ''
    title_var1_sc = ''
    # wet_models =['CIESM','CNRM-CM6-1-HR','CNRM-CM6-1','CNRM-ESM2-1','FGOALS-f3-L','INM-CM4-8','INM-CM5-0','NESM3'] #for annual non-scaled precipitation
    # dry_models =['ACCESS-ESM1-5','CESM2-WACCM','CESM2','CanESM5-CanOE','CanESM5','HadGEM3-GC31-LL','IPSL-CM6A-LR','KACE-1-0-G']
    # 30 and 70 percentiles
    wet_models = ['CIESM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','FGOALS-f3-L','INM-CM4-8','INM-CM5-0','MIROC6','NESM3']
    dry_models = ['ACCESS-ESM1-5','CESM2','CESM2-WACCM','CanESM5','CanESM5-CanOE','HadGEM3-GC31-LL','IPSL-CM6A-LR','KACE-1-0-G','MRI-ESM2-0','UKESM1-0-LL']


# Computations

# CMIP6

for season in season_list:

    # Select season:
    ds_period = time_period(ds,season)

    ds_rolled = ds_period.rolling(time=5,center=False).mean(dim='time') # 5 months rolling mean (Niño 3.4 definition)
    ds_MMmean = ds_rolled.mean(dim = 'model') # multimodel mean
    ds_MMmean_yearly = ds_MMmean.groupby(ds_MMmean.time.dt.year).mean(dim = 'time')
    ds_yearly = ds_rolled.groupby(ds_rolled.time.dt.year).mean(dim = 'time')
    wet_idx = np.zeros(len(wet_models))
    dry_idx = np.zeros(len(wet_models))

    for i in range(0,len(wet_models)):
        wet_idx[i] = np.where(ds_yearly.model_name == wet_models[i])[0]
        dry_idx[i] = np.where(ds_yearly.model_name == dry_models[i])[0]
    breakpoint()


    # COmpute 30 years trends
    years = np.linspace(2090,2100,11)
    m_cmip6_MMmean = []
    m_cmip6 = np.zeros((len(ds_yearly[:,0]),11))
    for i in range(0,11): #years
        m_cmip6_MMmean.append(np.polyfit(ds_MMmean_yearly.year[i:i+30],ds_MMmean_yearly[i:i+30], 1)[0])
       
        for j in range(0,len(ds_yearly[:,0])): #models
            m_cmip6[j,i] = np.polyfit(ds_yearly.year[i:i+30],ds_yearly[j,i:i+30], 1)[0]
    
    # Figure
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True)
    fig.suptitle(title_var1_sc+var+' 30 yr trend. '+title_sc)
    ax.set_title(title_ssp[ssp],loc='left',fontsize=9)
    ax.set_title(season_title[season],loc='right',fontsize=9)
    ax.scatter(m_cmip6[0,:],years, s = 5,color = 'blue',label = 'CMIP6 models mean')
    ax.scatter(m_cmip6[int(wet_idx[0]),:],years, s = 40,color = 'teal',label = 'Wet models')
    ax.scatter(m_cmip6[int(dry_idx[0]),:],years, s = 40,color = 'peru',label = 'Dry models')
    for i in range(0,len(m_cmip6[:,0])):
        ax.scatter(m_cmip6[i,:],years, s = 5,color = 'blue')
    for j in range(0,len(wet_idx)):
        i = wet_idx[j]
        ax.scatter(m_cmip6[int(i),:],years, s = 40,color = 'teal')#, marker = '>',label = ds_yearly.model_name[int(i)].values)#,color = 'teal')
        #ax.scatter(m_cmip6[int(i),:],years, s = 40,color = 'teal', marker = marker[j], label = ds_yearly.model_name[int(i)].values)#,color = 'teal')
    for j in range(0,len(dry_idx)):
        i = dry_idx[j]
        ax.scatter(m_cmip6[int(i),:],years, s = 40,color = 'peru')#, marker = 's',label = ds_yearly.model_name[int(i)].values)#,color = 'peru')
        #ax.scatter(m_cmip6[int(i),:],years, s = 40,color = 'peru', marker = marker[j],label = ds_yearly.model_name[int(i)].values)
    ax.scatter(m_cmip6_MMmean,years, s = 100,marker = 'd', color = 'navy',label = 'CMIP6 multimodel mean')
    ax.set_xlabel(clabel_trend[var])
    ax.set_ylabel('Year')
    ax.legend(loc = 'lower right')
    #fig.savefig('output_plots/'+var+sc+'_trend_future_CMIP6_observations.png')
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/'+var+sc+'_trend_CMIP6_observations_'+season+'_season_future.png')

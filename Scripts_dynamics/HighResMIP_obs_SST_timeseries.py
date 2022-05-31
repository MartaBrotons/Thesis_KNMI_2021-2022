# Python Script to compute 30 years trends for HighResMIP  models for the reference period (approx). (Figure A7)

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
from matplotlib.legend_handler import HandlerPatch


#------------------------INPUT------------------------#


trend = True

wdir = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/HighResMIP/03_preprocessed/HighResMIP/ts_Nino34/*ts*.nc'
# wdir = '/data/brotons/resampling/Method-2014/03_CMIP6_preprocessed/HighResMIP/ts_Nino34/'

colors = {'CMCC-CM2-HR4':'blue', 'CMCC-CM2-VHR4':'blue', 'CNRM-CM6-1':'purple', 'CNRM-CM6-1-HR':'purple', 'EC-Earth3P':'green', 'EC-Earth3P-HR':'green', 'HadGEM3-GC31-HH':'red', 'HadGEM3-GC31-HM':'red', 'HadGEM3-GC31-LL':'red', 'HadGEM3-GC31-MM':'red', 'MPI-ESM1-2-HR':'pink', 'MPI-ESM1-2-XR':'pink','CESM1-CAM5-SE-HR':'black', 'CESM1-CAM5-SE-LR':'black'}
markers = {'CMCC-CM2-HR4':'x', 'CMCC-CM2-VHR4':'d', 'CNRM-CM6-1':'x', 'CNRM-CM6-1-HR':'d', 'EC-Earth3P':'x', 'EC-Earth3P-HR':'d', 'HadGEM3-GC31-HH':'d', 'HadGEM3-GC31-HM':'o', 'HadGEM3-GC31-LL':'x', 'HadGEM3-GC31-MM':'o', 'MPI-ESM1-2-HR':'x', 'MPI-ESM1-2-XR':'d','CESM1-CAM5-SE-HR':'d', 'CESM1-CAM5-SE-LR':'x'}

# c = [ mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor=colors[i], edgecolor="none" ) for i in range(len(models))]
#patches = [ plt.plot([],[], marker=markers[i], ms=10, ls="", mec=None, color=colors[i], label=models[i]) for i in range(len(models))  ]
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

def boxmean_era(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.latitude))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('latitude','longitude'))#.dropna(dim='latitude').dropna(dim='longitude')
    return boxmean

def preprocessing_save_model_names(ds):
    ds = ds.assign_coords(model_name = ds.source_id)#.partition(' ')[0]) #to cut everything after the space in source
    return ds
# --------------------------------------------------------------------------#
breakpoint()

#HighResMIP data:
# datasets = []
# for filename in os.listdir(wdir):
#     if 'ts' in filename:
#         print(filename)
#         ds = xr.open_dataset(str(wdir)+str(filename))
#         datasets.append(ds)
# ts_HighResMIP = xr.concat(datasets,dim = 'model')  
breakpoint()
    
ts_HighResMIP = xr.open_mfdataset(wdir,concat_dim = 'model',combine = 'nested',preprocess = preprocessing_save_model_names )
ts_HighResMIP = ts_HighResMIP.load()
models = ts_HighResMIP.model_name.values
#ts_HighResMIP = ts_HighResMIP.sel(time = slice ('1-1-1979','31-12-2020'))
ts_HighResMIP_mean = boxmean(ts_HighResMIP.drop('time_bnds')) # if we don't drop the time_bnds variable sit gives a weird error 
ts_HighResMIP = ts_HighResMIP_mean.sel(time=slice('1-1-1979', '31-12-2020'),drop = True)
ts_HighResMIP_ref = ts_HighResMIP_mean.sel(time=slice('01-01-1991', '31-12-2020')).mean(dim = 'time')

nino34_HighResMIP = ts_HighResMIP-ts_HighResMIP_ref

nino34_HighResMIP_rolled = nino34_HighResMIP.rolling(time=5,center=False).mean(dim='time') # 5 months rolling mean (Niño 3.4 definition)
#nino34_HighResMIP_MMmean = nino34_HighResMIP_rolled.mean(dim = 'model') # multimodel mean
#nino34_HighResMIP_MMmean_yearly = nino34_HighResMIP_MMmean.groupby(ts_HighResMIP.time.dt.year).mean(dim = 'time')
nino34_HighResMIP_yearly = nino34_HighResMIP_rolled.groupby(ts_HighResMIP.time.dt.year).mean(dim = 'time')

breakpoint()

#ERA5 reanalysis
#sst_era5 = xr.open_dataset('/data/brotons/resampling/Method-2014/observations/ERA-5_SST_nino34.nc')
sst_era5 = xr.open_dataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_SST_glob.nc')
sst_era5 = sst_era5.where((sst_era5.latitude > -5) & (sst_era5.latitude < 5) & (sst_era5.longitude > 190) & (sst_era5.longitude < 240), drop=True)
sst_era5 = sst_era5.where(sst_era5.expver == 1,drop = True)
#sst_era5 = sst_era5.sel(time = slice ('1-1-1979','31-12-2020'))
sst_era5_mean = boxmean_era(sst_era5)
sst_era5 = sst_era5_mean.sel(time=slice('1-1-1979', '31-12-2020'),drop = True)
sst_era5_ref = sst_era5_mean.sel(time=slice('01-01-1991', '31-12-2020')).mean(dim = 'time')
nino34_era5 = sst_era5 -sst_era5_ref

nino34_era5_rolled = nino34_era5.rolling(time=5,center=False).mean(dim='time').dropna(dim='time')
nino34_era5_yearly = nino34_era5_rolled.groupby(nino34_era5_rolled.time.dt.year).mean(dim = 'time').dropna(dim='year')

breakpoint()
if trend == True:
    years = np.linspace(2010,2020,11)
    #m_HighResMIP_MMmean = []
    m_era5 = []
    m_HighResMIP = np.zeros((11,len(nino34_HighResMIP_yearly.model_name)))
    for i in range(0,11):
        #m_HighResMIP_MMmean.append(np.polyfit(nino34_HighResMIP_MMmean_yearly.year[i:i+30],nino34_HighResMIP_MMmean_yearly.ts[i:i+30], 1)[0])
        m_era5.append(np.polyfit(nino34_era5_yearly.year[i:i+30],nino34_era5_yearly.sst[i:i+30], 1)[0])
        for j in range(0,len(nino34_HighResMIP_yearly.model_name)):
            m_HighResMIP[i,j] = np.polyfit(nino34_HighResMIP_yearly.year[i:i+30],nino34_HighResMIP_yearly.ts[i:i+30,j], 1)[0]

    breakpoint()

    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(14,7),sharex=True, sharey =True)
    #ax.scatter(m_HighResMIP[:,0],years,color = colors[0],marker = markers[0],s=65,label = 'HighResMIP models mean')
    #legends = []
    for i in range(0,len(models)):
        model = models[i]
        ax.scatter(m_HighResMIP[:,i],years,s=65,color = colors[model],marker = markers[model],label=model)
        #legends.append(mlines.Line2D([], [], color=colors[i], label=models[i],marker = markers[i]))
    #ax.scatter(m_HighResMIP_MMmean,years, s = 65,color = 'navy',label = 'HighResMIP multimodel mean')
    ax.scatter(m_era5,years, s = 65,marker = 'X',color = 'orange',label = 'ERA5 reanalysis data')
    ax.set_xlabel('Niño 3.4 Trend [K/year]',fontsize = 14)
    ax.set_ylabel('Year',fontsize = 14)
    #ax.legend(c,models,bbox_to_anchor=(0.5, 0.5), loc='center', ncol=2, handler_map={mpatches.Circle: HandlerEllipse()}).get_frame().set_facecolor('#00FFCC')
    #ax.legend(bbox_to_anchor=(1.50, 0.7),loc='upper right', ncol=2, numpoints=1)#handles = [])
    ax.legend(fontsize = 14, loc='upper right', bbox_to_anchor=(1.3, 1.0),ncol=1, fancybox=True, shadow=True)
    fig.tight_layout()
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/Niño34_trend_HighResMIP_observations.png',dpi = 1200)


if trend == False: 
    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True)
    #ax.scatter(nino34_HighResMIP_yearly.ts[:,0],nino34_HighResMIP_yearly.year, s = 65,color = 'blue',label = 'HighResMIP models mean')
    for i in range(1,len(models)):
        model = models[i]
        ax.scatter(nino34_HighResMIP_yearly.ts[:,i],nino34_HighResMIP_yearly.year, s = 65,color = colors[model],marker = markers[model],label = model)
    #ax.scatter(nino34_HighResMIP_MMmean_yearly.ts,nino34_HighResMIP_MMmean_yearly.year, label = 'HighResMIP multimodel mean',color = 'navy',s=65)
    ax.scatter(nino34_era5_yearly.sst,nino34_era5_yearly.year,label = 'ERA5 reanalysis data',marker = 'X',color = 'orange',s=65)
    ax.set_xlabel('Niño 3.4 [K]')
    ax.set_ylabel('Year')
    ax.legend()
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/Niño34_HighResMIP_observations.png',dpi = 1200)

breakpoint()


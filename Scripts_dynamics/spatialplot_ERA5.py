# Script to compute 30 year climatologies for SST + wind (ERA5) and SLP (ERA5) + pr (GPCP). (Figures 3 and 4)
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
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import cmocean.cm as cmo
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

#------------------------INPUT------------------------#


# SETTINGS PROGRAM

month_list = [1,7] #[8,3]# [12,9]
month_title = ['January','July'] #['September','March'] ['December', 'September']

# coordinates domain
lat_matplot1 = [-20,70] 
lon_matplot1 = [50,250] 

lat_matplot2 = [-10,40] 
lon_matplot2 = [60,250] 

# -------------Things to choose:----------------------


# ----------------------------------------------------

paths = {'pr': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/gpcp_global.nc',
        'u': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_u10_glob.nc',
        'v': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_v10_glob.nc',
        'sst': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_SST_glob.nc',
        'slp': '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5/ERA5_mslp_glob.nc'}#} #'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/03_CMIP6_preprocessed/CMIP6/pr_NAtl_NPac/*ssp585*.nc' }#'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/ERA5_pr_global.nc'}


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
    weights = np.cos(np.deg2rad(da.latitude))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('latitude','longitude'))
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

#------------------------CALCULATIONS------------------------#
#breakpoint()
# Import variables
pr = xr.open_dataset(paths['pr'])#xr.open_mfdataset(paths['pr'], concat_dim = 'model', combine = 'nested', preprocess = preprocessing)
u10 = xr.open_dataset(paths['u'])
v10 = xr.open_dataset(paths['v'])
sst = xr.open_dataset(paths['sst'])
slp = xr.open_dataset(paths['slp'])


u10 = u10.mean(dim = 'expver',skipna = True)
v10 = v10.mean(dim = 'expver',skipna = True)
sst = sst.mean(dim = 'expver',skipna = True)
slp = slp.mean(dim = 'expver',skipna = True)

pr = pr.where((pr.time.dt.year>=1980) & (pr.time.dt.year<=2020),drop=True) 
u10 = u10.where((u10.time.dt.year>=1980) & (u10.time.dt.year<=2020),drop=True)
v10 = v10.where((v10.time.dt.year>=1980) & (v10.time.dt.year<=2020),drop=True)
sst = sst.where((sst.time.dt.year>=1980) & (sst.time.dt.year<=2020),drop=True)
slp = slp.where((slp.time.dt.year>=1980) & (slp.time.dt.year<=2020),drop=True)

# Start computations

for month,column in zip(month_list,[0,1]): 

    pr_season = pr.where((pr.time.dt.month>=month) & (pr.time.dt.month<=month),drop=True) #time_period(pr,season)
    pr_toplot = pr_season.mean(dim = 'time') #.mean(dim = 'model')

    slp_season = slp.where((slp.time.dt.month>=month) & (slp.time.dt.month<=month),drop=True) #time_period(slp,season)
    slp_toplot = slp_season.mean(dim = 'time')


    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    fig.suptitle('1980-2020 observations mean')#'ERA5 1982-2022 mean')
    ax.set_title(month_title[column],loc='right',fontsize=9)
    #ax.add_feature(cfeature.COASTLINE)
    #ax.coastlines()

    levels = np.arange(990,1040,4)
    z1_plot = ax.pcolormesh(pr_toplot.lon,pr_toplot.lat,pr_toplot.precip,transform=ccrs.PlateCarree(central_longitude=0), cmap=matplotlib.cm.get_cmap('cmo.rain' ))
    z1_plot.set_clim((0,10))
    cbar = fig.colorbar(z1_plot, ax=ax, shrink=0.8, orientation='horizontal')#'horizontal')
    cbar.set_label('Precipitation [mm/day]',fontsize = 14)
    cbar.ax.tick_params(labelsize = 14)
    z2_plot = ax.contour(slp_toplot.longitude,slp_toplot.latitude,slp_toplot.msl*1e-2,levels,transform=ccrs.PlateCarree(central_longitude=0),colors='k',linewidths=0.5)
    z2_plot.set_clim((9.5e2,1.5e3))
    ax.clabel(z2_plot, fontsize=6, inline=1,fmt= '% 1.1f')
    ax.coastlines(resolution = '50m')



    #ax.clabel(z2_plot, fontsize=6, inline=1,fmt= '% 1.1f')
    ax.set_xlim(lon_matplot1)
    ax.set_ylim(lat_matplot1)

    # # Define the xticks for longitude
    # ax.set_xticks(np.arange(lon_matplot1[0],lon_matplot1[1],20), crs=ccrs.PlateCarree())
    # lon_formatter = cticker.LongitudeFormatter()
    # ax.xaxis.set_major_formatter(lon_formatter)

    # # Define the yticks for latitude
    # ax.set_yticks(np.arange(lat_matplot1[0],lat_matplot1[1],20), crs=ccrs.PlateCarree())
    # lat_formatter = cticker.LatitudeFormatter()
    # ax.yaxis.set_major_formatter(lat_formatter)

    # fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/GPCP_ERA5_pr_slp_climatologies_'+month_title[column]+'.png')
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/GPCP_ERA5_pr_slp_climatologies_'+month_title[column]+'.png',dpi = 1200)
    #breakpoint()
#breakpoint()
for month,column in zip(month_list,[0,1]): 

    u10_season = u10.where((u10.time.dt.month>=month) & (u10.time.dt.month<=month),drop=True) #time_period(u10,season)
    u10_toplot = u10_season.mean(dim = 'time')
    v10_season = v10.where((v10.time.dt.month>=month) & (v10.time.dt.month<=month),drop=True) #time_period(v10,season)
    v10_toplot = v10_season.mean(dim = 'time')
    sst_season = sst.where((sst.time.dt.month>=month) & (sst.time.dt.month<=month),drop=True) #time_period(sst,season)
    sst_toplot = sst_season.mean(dim = 'time')

    fig , ax = plt.subplots(ncols=1, nrows=1,figsize=(10,7),sharex=True, sharey =True,subplot_kw=dict(projection = ccrs.PlateCarree(central_longitude=180)))
    fig.suptitle('1980-2020 observations mean')#'ERA5 1982-2022 mean')

    ax.set_title(month_title[column],loc='right',fontsize=9)
    # ax.add_feature(cfeature.COASTLINE)
    #ax.coastlines()
    #breakpoint()
    z3_plot = ax.pcolormesh(sst_toplot.longitude,sst_toplot.latitude,sst_toplot.sst-273.15,transform=ccrs.PlateCarree(central_longitude=0), cmap=matplotlib.cm.get_cmap('cmo.thermal' ))
    z4_plot = ax.contour(sst_toplot.longitude,sst_toplot.latitude,sst_toplot.sst-273.15,[0,28.5],transform=ccrs.PlateCarree(central_longitude=0),colors='white',linewidths=1.3)
    #z3_plot = ax.pcolormesh(u10_toplot.longitude, u10_toplot.latitude,(u10_toplot.u10**2+v10_toplot.v10**2)**(1/2),transform=ccrs.PlateCarree(central_longitude=0), cmap=matplotlib.cm.get_cmap('cmo.speed' ))
    z3_plot.set_clim((15,30))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='gray'))
    ax.coastlines(resolution = '50m')
    #z3_plot.set_clim((0,18))
    cbar = fig.colorbar(z3_plot, ax=ax, shrink=0.8, orientation='horizontal')
    cbar.set_label(r'SST [$^{\circ}C$]',fontsize = 14)
    cbar.ax.tick_params(labelsize = 14)
    z4_plot = ax.quiver(u10_toplot.longitude, u10_toplot.latitude, u10_toplot.u10.values, v10_toplot.v10.values,angles='xy', scale_units='xy', scale=1,width=0.0014,transform=ccrs.PlateCarree(central_longitude=0),regrid_shape=35)
    ax.quiverkey(z4_plot,0.04,-0.2,5, '5 m/s')
    ax.set_xlim(lon_matplot2)
    ax.set_ylim(lat_matplot2)

    # # Define the xticks for longitude
    # ax.set_xticks(np.arange(lon_matplot2[0],lon_matplot2[1],60), crs=ccrs.PlateCarree())
    # lon_formatter = cticker.LongitudeFormatter()
    # ax.xaxis.set_major_formatter(lon_formatter)

    # # Define the yticks for latitude
    # ax.set_yticks(np.arange(lat_matplot2[0],lat_matplot2[1],20), crs=ccrs.PlateCarree())
    # lat_formatter = cticker.LatitudeFormatter()
    # ax.yaxis.set_major_formatter(lat_formatter)

    
    # z1_plot.set_clim((0,10))
    # cbar = fig.colorbar(z1_plot, ax=ax, shrink=0.8, orientation='horizontal')
    # cbar.set_label('GPCP precipitation [mm/day]')
    #ax.coastlines(resolution = '50m')


    # fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/GPCP_ERA5_sst_wind_climatologies_'+month_title[column]+'.png')
    fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/GPCP_ERA5_sst_wind_climatologies_'+month_title[column]+'.png', dpi = 1200)
    #fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/GPCP_ERA5_wind_climatologies_'+month_title[column]+'.png')
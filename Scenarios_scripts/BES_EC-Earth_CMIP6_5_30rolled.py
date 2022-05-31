# This script plots temperature and precipitation timeseries from 16 EC-Earth output (5yr rolled) and 32 CMIP6 output (30yr rolled) for the wet and dry seasons (Figure 9)


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.lines as mlines
import glob
import netCDF4

#------------------------INPUT------------------------#

# Functions

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
    if 'lat' in ds:
        ds = ds.drop('lat')
    if 'lon' in ds:
        ds = ds.drop('lon')
    return ds.assign_coords(model_name = ds.source_id)

def boxmean(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.lat))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('lat','lon'))
    return boxmean

#------------------------CALCULATIONS------------------------#

# Import data (we already import lat lon mean data to run the code faster!)
tas_carib =  xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/ECE_bis/tas_BES/tas_BES_region_mean/tas_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble_member',combine = 'nested',preprocess = preprocessing)
tas_carib = tas_carib.load()
pr_carib =  xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/ECE_bis/pr_BES/pr_BES_region_mean/pr_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble_member',combine = 'nested',preprocess = preprocessing)
pr_carib = pr_carib.load()
pr_carib_cmip6 = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/CMIP6/pr_BES/pr_BES_region_mean/*ssp*5*.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing)
pr_carib_cmip6 = pr_carib_cmip6.load()
tas_carib_cmip6 = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/CMIP6/tas_BES/tas_BES_region_mean/*ssp*5*.nc',combine = 'nested', concat_dim = 'model',preprocess = preprocessing)
tas_carib_cmip6 = tas_carib_cmip6.load()

##EC-Earth:

#Temperature Caribbean
#Wet months:
tas_carib_wet = tas_carib.where((tas_carib.time.dt.month>=5) & (tas_carib.time.dt.month<=11)).dropna("time") #select months of the season
tas_carib_wet_ref = tas_carib_wet.tas.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') #average the temperature for the reference period (1991-2020)
tas_carib_wet = tas_carib_wet.rolling(time=5*7,center=False).mean(dim='time').dropna("time") #compute running average 5 yr window
tas_carib_wet_year = tas_carib_wet.groupby('time.year').mean(dim='time') #group by year
tas_carib_wet_year = tas_carib_wet_year.where((tas_carib_wet_year.year>=1991)&(tas_carib_wet_year.year<=2120)) #select period of time we are interested on 

#Dry months:
tas_carib_dry = tas_carib.where((tas_carib.time.dt.month<=4) | (tas_carib.time.dt.month>=12)).dropna("time")
tas_carib_dry_ref = tas_carib_dry.tas.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') 
tas_carib_dry = tas_carib_dry.rolling(time=5*5,center=False).mean(dim='time').dropna("time")
tas_carib_dry_year = tas_carib_dry.groupby('time.year').mean(dim='time')
tas_carib_dry_year = tas_carib_dry_year.where((tas_carib_dry_year.year>=1991)&(tas_carib_dry_year.year<=2120))

#Precipitation Caribbean
#Wet months:
pr_carib_wet = pr_carib.where((pr_carib.time.dt.month>=5)&(pr_carib.time.dt.month<=11)).dropna("time")
pr_carib_wet_ref = pr_carib_wet.pr.sel(time=slice('1-1-1991','31-12-2020')).mean(dim='time')
pr_carib_wet = pr_carib_wet.rolling(time=5*7,center=False).mean(dim='time').dropna("time")
pr_carib_wet_year = pr_carib_wet.groupby('time.year').mean(dim='time')
pr_carib_wet_year = pr_carib_wet_year.where((pr_carib_wet_year.year>=1991)&(pr_carib_wet_year.year<=2120))
#Dry months:
pr_carib_dry = pr_carib.where((pr_carib.time.dt.month<=4) | (pr_carib.time.dt.month>=12)).dropna("time")
pr_carib_dry_ref = pr_carib_dry.pr.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') #average the temperature for the reference period (1991-2020)
pr_carib_dry = pr_carib_dry.rolling(time=5*5,center=False).mean(dim='time').dropna("time")
pr_carib_dry_year = pr_carib_dry.groupby('time.year').mean(dim='time')
pr_carib_dry_year = pr_carib_dry_year.where((pr_carib_dry_year.year>=1991) & (pr_carib_dry_year.year<=2120) )

## CMIP6:

#Precipitation Caribbean
#Wet months:
pr_carib_cmip6_wet = pr_carib_cmip6.where((pr_carib_cmip6.time.dt.month>=5)&(pr_carib_cmip6.time.dt.month<=11)).dropna("time") 
pr_carib_cmip6_wet = pr_carib_cmip6_wet.rolling(time=30*7,center=False).mean(dim='time').dropna("time") #compute running average 30 yr window
pr_carib_cmip6_wet_year = pr_carib_cmip6_wet.groupby('time.year').mean(dim='time') 
pr_carib_cmip6_wet_ref = pr_carib_cmip6_wet.pr.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') 
pr_carib_cmip6_wet_year = pr_carib_cmip6_wet_year.where((pr_carib_cmip6_wet_year.year>=1991) & (pr_carib_cmip6_wet_year.year<=2120) )
#Dry months:
pr_carib_cmip6_dry = pr_carib_cmip6.where((pr_carib_cmip6.time.dt.month<=4) | (pr_carib_cmip6.time.dt.month>=12)).dropna("time")
pr_carib_cmip6_dry = pr_carib_cmip6_dry.rolling(time=30*5,center=False).mean(dim='time').dropna("time")
pr_carib_cmip6_dry_year = pr_carib_cmip6_dry.groupby('time.year').mean(dim='time')
pr_carib_cmip6_dry_ref = pr_carib_cmip6_dry.pr.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time')
pr_carib_cmip6_dry_year = pr_carib_cmip6_dry_year.where((pr_carib_cmip6_dry_year.year>=1991) & (pr_carib_cmip6_dry_year.year<=2120) )

#Temperature Caribbean
#Wet months:
tas_carib_cmip6_wet = tas_carib_cmip6.where((tas_carib_cmip6.time.dt.month>=5)&(tas_carib_cmip6.time.dt.month<=11)).dropna("time")
tas_carib_cmip6_wet = tas_carib_cmip6_wet.rolling(time=30*7,center=False).mean(dim='time').dropna("time") 
tas_carib_cmip6_wet_year = tas_carib_cmip6_wet.groupby('time.year').mean(dim='time')
tas_carib_cmip6_wet_ref = tas_carib_cmip6_wet.tas.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') 
tas_carib_cmip6_wet_year = tas_carib_cmip6_wet_year.where((tas_carib_cmip6_wet_year.year>=1991) & (tas_carib_cmip6_wet_year.year<=2120) )
#Dry months:
tas_carib_cmip6_dry = tas_carib_cmip6.where((tas_carib_cmip6.time.dt.month<=4) | (tas_carib_cmip6.time.dt.month>=12)).dropna("time")
tas_carib_cmip6_dry = tas_carib_cmip6_dry.rolling(time=30*5,center=False).mean(dim='time').dropna("time")
tas_carib_cmip6_dry_year = tas_carib_cmip6_dry.groupby('time.year').mean(dim='time')
tas_carib_cmip6_dry_ref = tas_carib_cmip6_dry.tas.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='time') 
tas_carib_cmip6_dry_year = tas_carib_cmip6_dry_year.where((tas_carib_cmip6_dry_year.year>=1991) & (tas_carib_cmip6_dry_year.year<=2120) )

#------------------------OUTPUT------------------------#

cmip6 = mlines.Line2D([], [], color='green', label='CMIP6') # to set manually the legend
ecearth = mlines.Line2D([], [], color='black', label='EC-Earth')

fig, ax = plt.subplots(2,2 ,figsize=(16,12))

for i in range(0,len(pr_carib_cmip6_dry_year.model)): # loop to plot one line timeseries per CMIP6 model
    ax[1,1].plot(tas_carib_cmip6_wet_year.year,tas_carib_cmip6_wet_year.tas[:,i,0,0]-tas_carib_cmip6_wet_ref[i,0,0], color='green')
    ax[1,0].plot(tas_carib_cmip6_dry_year.year,tas_carib_cmip6_dry_year.tas[:,i,0,0]-tas_carib_cmip6_dry_ref[i,0,0], color='green')
    ax[0,1].plot(pr_carib_cmip6_wet_year.year,(pr_carib_cmip6_wet_year.pr[:,i,0,0]-pr_carib_cmip6_wet_ref[i,0,0])*(100/pr_carib_cmip6_wet_ref[i,0,0]), color='green')
    ax[0,0].plot(pr_carib_cmip6_dry_year.year,(pr_carib_cmip6_dry_year.pr[:,i,0,0]-pr_carib_cmip6_dry_ref[i,0,0])*(100/pr_carib_cmip6_dry_ref[i,0,0]), color='green')

for i in range(0,len(pr_carib_dry_year.ensemble_member)): # loop to plot one line timeseries per EC-Earth ensemble
    ax[1,1].plot(tas_carib_wet_year.year,tas_carib_wet_year.tas[:,i,0,0]-tas_carib_wet_ref[i,0,0], color='black')
    ax[1,0].plot(tas_carib_dry_year.year,tas_carib_dry_year.tas[:,i,0,0]-tas_carib_dry_ref[i,0,0], color='black')
    ax[0,1].plot(pr_carib_wet_year.year,(pr_carib_wet_year.pr[:,i,0,0]-pr_carib_wet_ref[i,0,0])*(100/pr_carib_wet_ref[i,0,0]),color='black')
    ax[0,0].plot(pr_carib_dry_year.year,(pr_carib_dry_year.pr[:,i,0,0]-pr_carib_dry_ref[i,0,0])*(100/pr_carib_dry_ref[i,0,0]),color='black')

ax[1,1].grid(True)
ax[1,1].set_ylim(-1,5.2)
ax[1,1].set_xlim(1990,2120)
ax[1,1].set_xlabel('Year',fontsize = 16)
ax[1,1].set_ylabel(r'$\Delta$ Temperature  [$^{\circ} C$]',fontsize = 16)
ax[1,1].set_title('Wet season (May-November)',loc = 'right',fontsize = 14) 
plt.xticks(fontsize=12)

ax[1,0].grid(True)
ax[1,0].set_ylim(-1,5.2)
ax[1,0].set_xlim(1990,2120)
ax[1,0].set_xlabel('Year',fontsize = 16)
ax[1,0].set_ylabel(r'$\Delta$ Temperature [$^{\circ} C$]',fontsize = 16)
ax[1,0].set_title('Dry season (December-April)',loc = 'right',fontsize = 14) 
plt.xticks(fontsize=12)


ax[0,1].grid(True)
ax[0,1].set_xlabel('Year',fontsize = 16)
ax[0,1].set_xlim(1990,2120)
ax[0,1].set_ylabel(r'$\Delta$ Precipitation $[\%]$',fontsize = 16)
ax[0,1].set_title('Wet season (May-November)',loc = 'right',fontsize = 14) 
ax[0,1].set_ylim(-50,30)


ax[0,0].grid(True)
ax[0,0].set_xlabel('Year',fontsize = 16)
ax[0,0].set_xlim(1990,2120)
ax[0,0].set_ylim(-50,30)
ax[0,0].set_ylabel(r'$\Delta$ Precipitation $[\%]$',fontsize = 16)
ax[0,0].set_title('Dry season (December-April)',loc = 'right',fontsize = 14) 
ax[1,1].legend(handles=[cmip6,ecearth],fontsize = 16)

ax[1,1].tick_params(axis='both', which='major', labelsize=14) # to set all the x and y ticks to the same size
ax[1,0].tick_params(axis='both', which='major', labelsize=14)
ax[0,1].tick_params(axis='both', which='major', labelsize=14)
ax[0,0].tick_params(axis='both', which='major', labelsize=14)


fig.suptitle('EC-Earth (5 yr running mean) + CMIP6 (30 yr running mean)',fontsize=20)
fig.savefig('/data/brotons/resampling/Method-2014/output_plots/BES_EC-Earth_CMIP6_5_30rolled.png')
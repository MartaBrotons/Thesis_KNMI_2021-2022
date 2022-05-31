##!/usr/bin/env python3
#Python script to compute and plot CMIP-6 bandbreedte in tas, pr and ev relative or absolute change
#can be run with e.g. conda create -n py_envir xarray netCDF4 ipython pylint matplotlib

#Script to compute the statistics without lat lon averages for the tropics region. Output is saved on a folder

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
from cftime import DatetimeNoLeap
from tqdm import tqdm

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


# CMIP6
breakpoint()

var_list = ['tas','psl','hfls','ts','ua200','ua850','va200','va850','wap850','zg850','moist_adv','moist_conv','moist_mass_conv'] #['pr', 'tas','psl','hfls','hfss','hur850','hur925','hus850','hus925','ta200','ta850','ts','ua200','ua850','va200','va850','wap500','wap850','zg850']

for var in var_list:

    if (var == 'pr') or (var == 'tas') or (var == 'psl'): starting_date = '1850-01-16'
    else: starting_date = '1950-01-16'

    wdir_pr = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/'+var+'/'
    outdir = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/02_time_selected/'+var+'/'

    #to set all the cmip6/ec-earth models/runs to the same date time coordinate:

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for filename in tqdm(os.listdir(wdir_pr)): #for every file in the directory
        ds=xr.open_dataset(str(wdir_pr)+str(filename), decode_times = True) #open .nc file
        if ds.indexes['time'].inferred_type != "datetime64":
            ds['time'] = ds.indexes['time'].to_datetimeindex()
        ds['time'] = pd.date_range(start=starting_date,end='2100-12-16',periods=len(ds.time))
        print(filename)

        ds = ds.assign_coords(model_name = filename.partition('_Amon')[0])
        ds_dry = time_period(ds,'dry')
        ds_2005_dry = ds_dry.where(((ds_dry.time.dt.year >= 1991 ) & (ds_dry.time.dt.year <= 2020)), drop=True).mean(dim = 'time')
        ds_2085_dry = ds_dry.where(((ds_dry.time.dt.year >= 2071 ) & (ds_dry.time.dt.year <= 2100)) , drop=True).mean(dim = 'time')
        ds_wet = time_period(ds,'wet')
        ds_2005_wet = ds_wet.where(((ds_wet.time.dt.year >= 1991 ) & (ds_wet.time.dt.year <= 2020)), drop=True).mean(dim = 'time')
        ds_2085_wet = ds_wet.where(((ds_wet.time.dt.year >= 2071 ) & (ds_wet.time.dt.year <= 2100)) , drop=True).mean(dim = 'time')

        ds_2005_dry.to_netcdf(outdir+filename.partition('.nc')[0]+'_dry_season_2005.nc')
        ds_2085_dry.to_netcdf(outdir+filename.partition('.nc')[0]+'_dry_season_2085.nc')
        ds_2005_wet.to_netcdf(outdir+filename.partition('.nc')[0]+'_wet_season_2005.nc')
        ds_2085_wet.to_netcdf(outdir+filename.partition('.nc')[0]+'_wet_season_2085.nc')



# HighResMIP

# wdir_pr = '/data/brotons/resampling/Method-2014/00_join_CMIP6_models/HighResMIP/ts/'
# outdir = '/data/brotons/resampling/Method-2014/01_datetime_CMIP6_models/HighResMIP/ts/'

# #to set all the cmip6/ec-earth models/runs to the same date time coordinate:

# if not os.path.exists(outdir):
#     os.makedirs(outdir)


# for filename in tqdm(os.listdir(wdir_pr)): #for every file in the directory
#     ds=xr.open_dataset(str(wdir_pr)+str(filename), decode_times = True) #open .nc file
#     if ds.indexes['time'].inferred_type != "datetime64":
#         ds['time'] = ds.indexes['time'].to_datetimeindex()
#     ds['time'] = pd.date_range(start='1950-01-16',end='2050-12-16',periods=len(ds.time))
#     print(filename)
#     ds.to_netcdf(outdir+filename)



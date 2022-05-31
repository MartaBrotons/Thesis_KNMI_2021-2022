# Script to compute the scenarios mean for annual, dry and wet seasons

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

# ------------------ Things to select ----------------------

target_year = '2050'

#list of things to do:

scenarios_list = ['GL', 'GH', 'WL', 'WH']
season_list = ['annual', 'dry', 'wet']

# ----------------------------------------------------------


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


#------------------------INPUT------------------------#

path = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/Scenarios_output/resampled_' + target_year + '/'

dtas = np.empty((len(season_list),len(scenarios_list)))
dtas[:,:] = np.nan
dpr = np.empty((len(season_list),len(scenarios_list)))
dpr[:,:] = np.nan

tas_control = np.empty((len(season_list),len(scenarios_list)))
tas_control[:,:] = np.nan
pr_control = np.empty((len(season_list),len(scenarios_list)))
pr_control[:,:] = np.nan

for scenario,j in zip(scenarios_list,range(0,len(scenarios_list))):
    control = xr.open_dataset(path + 'resampled_control_' + scenario + '_EOC.nc')
    future = xr.open_dataset(path + 'resampled_future_' + scenario + '_EOC.nc')

    for season,i in zip(season_list,range(0,len(season_list))): 

        control_period = time_period(control, season)
        future_period = time_period(future, season)

        tas_control[i,j] = control_period.tas.mean(dim = 'time').mean(dim = 'sample')-273.15
        pr_control[i,j] = control_period.pr.mean(dim = 'time').mean(dim = 'sample')* 86400

        dpr[i,j] = ((future_period.pr.mean(dim = 'time').mean(dim = 'sample') - control_period.pr.mean(dim = 'time').mean(dim = 'sample'))/control_period.pr.mean(dim = 'time').mean(dim = 'sample'))*100
        dtas[i,j] = future_period.tas.mean(dim = 'time').mean(dim = 'sample') - control_period.tas.mean(dim = 'time').mean(dim = 'sample')

print('tas control')
print(tas_control )

print('pr control')
print(pr_control )

print('dtas')
print(dtas )

print('dpr')
print(dpr )

pd.DataFrame(tas_control).to_csv(path+ 'tas_control.csv', header = scenarios_list)
pd.DataFrame(pr_control).to_csv(path+ 'pr_control.csv', header = scenarios_list)

pd.DataFrame(dtas).to_csv(path+ 'dtas.csv', header = scenarios_list)
pd.DataFrame(dpr).to_csv(path+ 'dpr.csv', header = scenarios_list)
##!/usr/bin/env python3
#Python script to compute and save the moisture massconvergence as in Martinez 2019 figure 2e
# Computed as 1000-7000 hPa. High topography regions (no 100hPa values) are masked.

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
from cftime import DatetimeNoLeap
from tqdm import tqdm



def divergence(u, v, lon, lat):
    '''calculate divergence'''
    cos = np.cos
    pi = np.pi
    a_earth    = 6371e3 # radius earth
    dudx       = (u.differentiate('lon')) * (180/pi)*((a_earth*cos(lat*pi/180))**(-1))
    vcoslat    = (v*cos(lat*pi/180))
    dvdy       = (vcoslat.differentiate('lat'))*(180/pi) * ((a_earth*cos(lat*pi/180))**(-1))
    divergence = (dudx + dvdy)
    divergence.attrs['long_name'] = 'divergence'
    divergence.attrs['units']     = 's**-1'
    divergence.attrs['history']   = 'computed as dudx + dvdy from u and v'
    return divergence


# CMIP6
breakpoint()


path = {'hus700':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/hus700/',
        'ua700':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/ua700/',
        'va700':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/va700/',
        'hus850':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/hus850/',
        'ua850':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/ua850/',
        'va850':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/va850/',
        'hus925':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/hus925/',
        'ua925':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/ua925/',
        'va925':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/va925/',
        'hus1000':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/hus1000/',
        'ua1000':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/ua1000/',
        'va1000':'/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/va1000/'}

    #to set all the cmip6/ec-earth models/runs to the same date time coordinate:

for model in tqdm(['EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorESM2-MM','UKESM1-0-LL']):#'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CanESM5','CanESM5-CanOE','CESM2','CESM2-WACCM','CIESM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorESM2-MM','UKESM1-0-LL']):
    ua1000 = xr.open_dataset(path['ua1000']+model+'_Amon_historical-ssp585_195001-210012_ua1000.nc')
    va1000 = xr.open_dataset(path['va1000']+model+'_Amon_historical-ssp585_195001-210012_va1000.nc')
    q1000 = xr.open_dataset(path['hus1000']+model+'_Amon_historical-ssp585_195001-210012_hus1000.nc')
    div_wind1000 = divergence(ua1000.ua,va1000.va,ua1000.lon,ua1000.lat)
    q_div_wind1000 = q1000.hus*div_wind1000
    ua925 = xr.open_dataset(path['ua925']+model+'_Amon_historical-ssp585_195001-210012_ua925.nc')
    va925 = xr.open_dataset(path['va925']+model+'_Amon_historical-ssp585_195001-210012_va925.nc')
    q925 = xr.open_dataset(path['hus925']+model+'_Amon_historical-ssp585_195001-210012_hus925.nc')
    div_wind925 = divergence(ua925.ua,va925.va,ua925.lon,ua925.lat)
    q_div_wind925 = q925.hus*div_wind925
    ua850 = xr.open_dataset(path['ua850']+model+'_Amon_historical-ssp585_195001-210012_ua850.nc')
    va850 = xr.open_dataset(path['va850']+model+'_Amon_historical-ssp585_195001-210012_va850.nc')
    q850 = xr.open_dataset(path['hus850']+model+'_Amon_historical-ssp585_195001-210012_hus850.nc')
    div_wind850 = divergence(ua850.ua,va850.va,ua850.lon,ua850.lat)
    q_div_wind850 = q850.hus*div_wind850
    ua700 = xr.open_dataset(path['ua700']+model+'_Amon_historical-ssp585_195001-210012_ua700.nc')
    va700 = xr.open_dataset(path['va700']+model+'_Amon_historical-ssp585_195001-210012_va700.nc')
    q700 = xr.open_dataset(path['hus700']+model+'_Amon_historical-ssp585_195001-210012_hus700.nc')
    div_wind700 = divergence(ua700.ua,va700.va,ua700.lon,ua700.lat)
    q_div_wind700 = q700.hus*div_wind700

    mass_conv = -(1/(9.8*1e3))*(((q_div_wind1000.drop('plev')+q_div_wind925.drop('plev'))/2)*(1000e2-925e2)+((q_div_wind925.drop('plev')+q_div_wind850.drop('plev'))/2)*(925e2-850e2)+((q_div_wind850.drop('plev')+q_div_wind700.drop('plev'))/2)*(850e2-700e2)) # in m/s
    mass_conv = mass_conv*1e3*86400 #trasnform to mm/day

    mass_conv = mass_conv.to_dataset(name = 'mass_conv')

    mass_conv['model_name'] = ua1000.source_id
    mass_conv['variant_label'] = ua1000.variant_label
    mass_conv['units'] = 'mm/day'
    mass_conv['long_name'] = 'moisture mass convergence'
    mass_conv.to_netcdf('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/moist_mass_conv/'+model+'_Amon_historical-ssp585_195001-210012_moist_mass_conv.nc')
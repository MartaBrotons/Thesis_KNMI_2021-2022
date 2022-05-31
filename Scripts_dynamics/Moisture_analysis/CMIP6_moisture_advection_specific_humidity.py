##!/usr/bin/env python3
#Python script to compute and save the advection of specific moisture by the mean flow as in Martinez 2019 figure 2f
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

def gradient(phi, lon, lat):
    '''calculate divergence'''
    cos = np.cos
    pi = np.pi
    a_earth    = 6371e3 # radius earth
    dphidx       = (phi.differentiate('lon')) * (180/pi)*((a_earth*cos(lat*pi/180))**(-1))
    phicoslat    = (phi*cos(lat*pi/180))
    dphidy       = (phicoslat.differentiate('lat'))*(180/pi) * ((a_earth*cos(lat*pi/180))**(-1))
    grad_phix = dphidx
    grad_phiy = dphidy
    grad_phix.attrs['long_name'] = 'Gradient x axis'; grad_phiy.attrs['long_name'] = 'Gradient y axis'
    grad_phix.attrs['units']     = 's**-1'; grad_phiy.attrs['units']     = 's**-1'
    return grad_phix, grad_phiy 


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

for model in tqdm(['ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CanESM5','CanESM5-CanOE','CESM2','CESM2-WACCM','CIESM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorESM2-MM','UKESM1-0-LL']):
    ua1000 = xr.open_dataset(path['ua1000']+model+'_Amon_historical-ssp585_195001-210012_ua1000.nc')
    va1000 = xr.open_dataset(path['va1000']+model+'_Amon_historical-ssp585_195001-210012_va1000.nc')
    q1000 = xr.open_dataset(path['hus1000']+model+'_Amon_historical-ssp585_195001-210012_hus1000.nc')
    gradx_q1000,grady_q1000 = gradient(q1000.hus,q1000.lon,q1000.lat)
    gradx_q1000 = gradx_q1000.drop('plev'); grady_q1000 = grady_q1000.drop('plev')
    ua_coords = ua1000.ua[:,:,:,:]
    adv1000 = np.empty((len(ua1000.time),len(ua1000.plev),len(ua1000.lat),len(ua1000.lon)))
    adv1000[:,:,:,:] = np.nan
    adv1000 = xr.DataArray(adv1000, coords=ua_coords.coords)
    ua925 = xr.open_dataset(path['ua925']+model+'_Amon_historical-ssp585_195001-210012_ua925.nc')
    va925 = xr.open_dataset(path['va925']+model+'_Amon_historical-ssp585_195001-210012_va925.nc')
    q925 = xr.open_dataset(path['hus925']+model+'_Amon_historical-ssp585_195001-210012_hus925.nc')
    div925 = divergence(q925.hus*ua925.ua,q925.hus*va925.va,ua925.lon,ua925.lat)
    gradx_q925,grady_q925 = gradient(q925.hus,q925.lon,q925.lat)
    gradx_q925 = gradx_q925.drop('plev'); grady_q925 = grady_q925.drop('plev')
    adv925 = np.empty((len(ua925.time),len(ua925.plev),len(ua925.lat),len(ua925.lon)))
    adv925[:,:,:,:] = np.nan
    adv925 = xr.DataArray(adv925, coords=ua_coords.coords)
    ua850 = xr.open_dataset(path['ua850']+model+'_Amon_historical-ssp585_195001-210012_ua850.nc')
    va850 = xr.open_dataset(path['va850']+model+'_Amon_historical-ssp585_195001-210012_va850.nc')
    q850 = xr.open_dataset(path['hus850']+model+'_Amon_historical-ssp585_195001-210012_hus850.nc')
    gradx_q850,grady_q850 = gradient(q850.hus,q850.lon,q850.lat)
    gradx_q850 = gradx_q850.drop('plev'); grady_q850 = grady_q850.drop('plev')
    adv850 = np.empty((len(ua850.time),len(ua850.plev),len(ua850.lat),len(ua850.lon)))
    adv850[:,:,:,:] = np.nan
    adv850 = xr.DataArray(adv850, coords=ua_coords.coords)
    ua700 = xr.open_dataset(path['ua700']+model+'_Amon_historical-ssp585_195001-210012_ua700.nc')
    va700 = xr.open_dataset(path['va700']+model+'_Amon_historical-ssp585_195001-210012_va700.nc')
    q700 = xr.open_dataset(path['hus700']+model+'_Amon_historical-ssp585_195001-210012_hus700.nc')
    gradx_q700,grady_q700 = gradient(q700.hus,q700.lon,q700.lat)
    gradx_q700 = gradx_q700.drop('plev'); grady_q700 = grady_q700.drop('plev')
    adv700 = np.empty((len(ua700.time),len(ua700.plev),len(ua700.lat),len(ua700.lon)))
    adv700[:,:,:,:] = np.nan
    adv700 = xr.DataArray(adv700, coords=ua_coords.coords)
    adv1000 = ua1000.ua*gradx_q1000 + va1000.va*gradx_q1000
    adv925 = ua925.ua*gradx_q925 + va925.va*gradx_q925
    adv850 = ua850.ua*gradx_q850 + va850.va*gradx_q850
    adv700 = ua700.ua*gradx_q700 + va700.va*gradx_q700

    adv = -(1/(9.8*1e3))*(((adv1000.drop('plev')+adv925.drop('plev'))/2)*(1000e2-925e2)+((adv925.drop('plev')+adv850.drop('plev'))/2)*(925e2-850e2)+((adv850.drop('plev')+adv700.drop('plev'))/2)*(850e2-700e2)) # in m/s
    adv = adv*1e3*86400 #transform to mm/day
    adv = adv.to_dataset(name = 'adv')

    adv['model_name'] = ua1000.source_id
    adv['variant_label'] = ua1000.variant_label
    adv['units'] = 'mm/day'
    adv['long_name'] = 'advection of specific moisture'
    adv.to_netcdf('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/moist_adv/'+model+'_Amon_historical-ssp585_195001-210012_moist_adv.nc')
#!/usr/bin/env python3
# Script to select a region from the whole globe, take the lat lon average, and save it into sepparated netcdf files (1 per model and forcing scenario)
# if input is ts on boxmean(ds.ts)

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
from cftime import DatetimeNoLeap
from tqdm import tqdm

def boxmean(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.lat))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('lat','lon'))
    return boxmean

def preprocessing(ds):
    ds = ds.assign_coords(model_name = ds.source_id)
    ds = ds.assign_coords(variant_label = ds.variant_label)
    return ds

domain = 'BES' # 'NAtl_NPac' '40N_40S'
var_list = ['tas'] #['pr','tas'] #['pr', 'tas','psl','hfls','hfss','hur850','hur925','hus850','hus925','ta200','ta850','ts','ua200','ua850','va200','va850','wap500','wap850','zg850']

breakpoint()
latit = {'BES': [5,20], 'NAtl_NPac': [-20,60],'40N_40S': [-40,40],'Nino34': [-5,5],'Nino12':[0,10],'NASPG':[50,60],'tropics':[-10,40],'NH_poles_nonNASPG':[50,60],'BES2': [10,20]}
longit = {'BES': [270,310], 'NAtl_NPac': [110,360],'40N_40S': [0,360],'Nino34':[120,170],'Nino12':[270,280],'NASPG':[305,340],'tropics':[230,340],'NH_poles_nonNASPG':[340,305],'BES2': [275,300]}

latitude = latit[domain]
longitude = longit[domain]

for var in var_list:
    indir = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/01_regrid_CMIP6_models/'+var+'/'
    outdir = '/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/'+domain+'/'+var+'/'
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    for filename in tqdm(os.listdir(indir)):
        ds = xr.open_mfdataset(indir+filename,preprocess = preprocessing)
        if domain == '40N_40S': ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]), drop=True)
        else: ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]) & (ds.lon > longitude[0]) & (ds.lon < longitude[1]), drop=True)
        ds = boxmean(ds) #ds.ts
        ds.to_netcdf(outdir+filename)

        
    # for filename in tqdm(os.listdir(indir)):
    #     ds = xr.open_mfdataset(indir+filename,preprocess = preprocessing)
    #     if domain == '40N_40S': ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]), drop=True)
    #     else:
    #         if longitude[0] <= longitude[1]: ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]) & (ds.lon > longitude[0]) & (ds.lon < longitude[1]), drop=True)
    #         elif longitude[0] >= longitude[1]:
    #             ds = ds.where((ds.lat > latitude[0]) & (ds.lat < latitude[1]),drop = True)
    #             ds = ds.where((ds.lon > longitude[0]) | (ds.lon < longitude[1]), drop=True)
    #     ds = boxmean(ds)
    #     ds.to_netcdf(outdir+filename)

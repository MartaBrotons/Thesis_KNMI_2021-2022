# This script generates rainfall climatologies on the region of the caribbean from observations (GPCP) and CMIP6 multimodel mean and 90% spread and EC-Earth3. (Figure 5 report)

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset
import pandas as pd



#function
def boxmean(da):
    """Compute spatial mean"""
    weights = np.cos(np.deg2rad(da.lat))
    weights.name = 'weights'
    boxmean = da.weighted(weights).mean(dim=('lat','lon'))
    return boxmean

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

#  Import data : 

#gpcp observations
gpcp = xr.open_dataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/observations/gpcp_global.nc')
gpcp_year = gpcp.sel(time=slice('1-1-1991', '31-12-2020'))
gpcp_month = gpcp_year.groupby(gpcp_year.time.dt.month).mean(dim = 'time')
gpcp_mean = boxmean(gpcp_month.where((gpcp_month.lat > 5) & (gpcp_month.lat < 20) & (gpcp_month.lon > 270) & (gpcp_month.lon < 310)))

breakpoint()
#ec-earth
ec_earth = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/03_CMIP6_preprocessed/ECE_bis/pr_BES/pr_BES_region_mean/pr*', concat_dim = 'ensemble', combine = 'nested')
ec_earth = ec_earth.load()
ec_earth_year = ec_earth.sel(time=slice('1-1-1991', '31-12-2020'))
ec_earth_month = ec_earth_year.groupby(ec_earth_year.time.dt.month).mean(dim = 'time')
ec_earth_mean = ec_earth_month.mean(dim = 'ensemble')

#cmip6
cmip6 = xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/CMIP6/03_area_mean/BES/pr/*', concat_dim = 'model', combine = 'nested', preprocess = preprocessing_area_mean)
cmip6 = cmip6.load()
cmip6_year = cmip6.sel(time=slice('1-1-1991', '31-12-2020'))
cmip6_month = cmip6_year.groupby(cmip6_year.time.dt.month).mean(dim = 'time')
cmip6_mean = cmip6_month.mean(dim = 'model')



months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec']

fig,ax = plt.subplots(1,1,figsize=(14,7))
plt.plot(gpcp_mean.month, gpcp_mean.precip, label='Observations GPCP', color = 'black',linewidth=3)
plt.plot(ec_earth_mean.month, ec_earth_mean.pr[:,0,0]*86400, label='EC-Earth3', color = 'red',linewidth=2)
plt.plot(cmip6_mean.month, cmip6_mean.pr*86400, label='CMIP6', color = 'blue',linewidth=2)

plt.legend(fontsize = 18,loc = 'upper left')
ax.vlines(4.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
ax.vlines(6.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
ax.vlines(9.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
ax.vlines(11.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
ax.text(2,0.3,'WDS',fontsize = 18,color = 'gray')
ax.text(5.2,0.3,'ERS',fontsize = 18,color = 'gray')
ax.text(7.75,0.3,'MSD',fontsize = 18,color = 'gray')
ax.text(10.2,0.3,'LRS',fontsize = 18,color = 'gray')
# ax.vlines(4.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
# ax.vlines(11.5,-1,8, color = 'gray',linestyle = '--', linewidth = 1)
# ax.text(1.5,0.3,'Dry season',fontsize = 18,color = 'gray')
# ax.text(7,0.3,'Wet season',fontsize = 18,color = 'gray')
ax.set_xticks(gpcp_mean.month)
ax.set_xticklabels(months,rotation=50,fontsize = 18)
ax.set_yticklabels(np.linspace(0,6,7),fontsize = 18)
ax.set_ylim(0,6.5)
#plt.xlabel("Month of the year")
plt.ylabel("Precipitation [mm/day]",fontsize = 18)
#plt.title("Caribbean precipitation climatologies (1991-2020)",fontsize = 20)
plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/BES_rainfall_climatologies_4seasons.png',dpi = 300)
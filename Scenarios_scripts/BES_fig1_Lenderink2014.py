# This script tries to generate figure 1 from Lenderink et al., 2014 for the Caribbean region. As input it need EC-Earth data
# Here we take the rolling mean BEFORE spliting our dataset into seasons (Figure 8)

#------------------------INPUT------------------------#

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


#------------------------CALCULATIONS------------------------#

# Import data (we import already lat lon average data to run the script faster!)
tas_carib =  xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/ECE_bis/tas_BES/tas_BES_region_mean/tas_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble_member',combine = 'nested')
tas_carib = tas_carib.load()
pr_carib =  xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/03_preprocessed/ECE_bis/pr_BES/pr_BES_region_mean/pr_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble_member',combine = 'nested')
pr_carib = pr_carib.load()
tas_glob =  xr.open_mfdataset('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/preprocessed_data/ECEbis/04_global_mean/ECEbis/tas/ssp585/tas_mon_EC-Earth3bis_historical-ssp585*.nc',concat_dim = 'ensemble_member',combine = 'nested')
tas_glob = tas_glob.load()

breakpoint()

tas_glob = tas_glob.tas[:,:,0,0] #skip lat and lon dimensions
tas_carib = tas_carib.tas[:,:,0,0]
pr_carib = pr_carib.pr[:,:,0,0]

tas_carib= tas_carib.sel(time=slice('1-1-1950', '31-12-2110')) # Select the time we want
tas_glob= tas_glob.sel(time=slice('1-1-1950', '31-12-2110'))
pr_carib= pr_carib.sel(time=slice('1-1-1950', '31-12-2110'))



#TEMPERATURE CARIBBEAN

#Wet months:
tas_carib_wet = tas_carib.where((tas_carib.time.dt.month>=5) & (tas_carib.time.dt.month<=11)).dropna("time") # select months of the season
tas_carib_wet_ref_mean = tas_carib_wet.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='ensemble_member').mean(dim='time') #average the temperature for the reference period (1991-2020)
tas_carib_wet = tas_carib_wet.rolling(time=30*7,center=False).mean(dim='time').dropna("time") #compute running average 30 yr window
tas_carib_wet_mean = tas_carib_wet.mean(dim='ensemble_member') #compute the mean for the wet season
tas_carib_wet_quant = tas_carib_wet.quantile([0.05,0.95],dim='ensemble_member') #compute the percentiles for the we season

#Dry months:
tas_carib_dry = tas_carib.where((tas_carib.time.dt.month<=4) | (tas_carib.time.dt.month>=12)).dropna("time")
tas_carib_dry_ref_mean = tas_carib_dry.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='ensemble_member').mean(dim='time') 
tas_carib_dry = tas_carib_dry.rolling(time=30*5,center=False).mean(dim='time').dropna("time")
tas_carib_dry_mean = tas_carib_dry.mean(dim='ensemble_member') 
tas_carib_dry_quant = tas_carib_dry.quantile([0.05,0.95],dim='ensemble_member')

#TEMPERATURE GLOBAL

#Wet months:
tas_glob_wet = tas_glob.where((tas_glob.time.dt.month>=5)&(tas_glob.time.dt.month<=11)).dropna("time")
tas_glob_wet_ref_mean = tas_glob_wet.sel(time=slice('1-1-1991','31-12-2020')).mean(dim='ensemble_member').mean(dim='time')
tas_glob_wet = tas_glob_wet.rolling(time=30*7,center=False).mean(dim='time').dropna("time")
tas_glob_wet_mean = tas_glob_wet.mean(dim='ensemble_member')
tas_glob_wet_quant = tas_glob_wet.quantile([0.05,0.95],dim='ensemble_member')

#Dry months:
tas_glob_dry = tas_glob.where((tas_glob.time.dt.month<=4) | (tas_glob.time.dt.month>=12)).dropna("time")
tas_glob_dry_ref_mean = tas_glob_dry.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='ensemble_member').mean(dim='time')
tas_glob_dry = tas_glob_dry.rolling(time=30*5,center=False).mean(dim='time').dropna("time")
tas_glob_dry_mean = tas_glob_dry.mean(dim='ensemble_member') 
tas_glob_dry_quant = tas_glob_dry.quantile([0.05,0.95],dim='ensemble_member')

#PRECIPITATION CARIBBEAN

#Wet months:
pr_carib_wet = pr_carib.where((pr_carib.time.dt.month>=5)&(pr_carib.time.dt.month<=11)).dropna("time")
pr_carib_wet_ref_mean = pr_carib_wet.sel(time=slice('1-1-1991','31-12-2020')).mean(dim='ensemble_member').mean(dim='time')
pr_carib_wet = pr_carib_wet.rolling(time=30*7,center=False).mean(dim='time').dropna("time")
pr_carib_wet_mean = pr_carib_wet.mean(dim='ensemble_member')
pr_carib_wet_quant = pr_carib_wet.quantile([0.05,0.95],dim='ensemble_member')

#Dry months:
pr_carib_dry = pr_carib.where((pr_carib.time.dt.month<=4) | (pr_carib.time.dt.month>=12)).dropna("time")
pr_carib_dry_ref_mean = pr_carib_dry.sel(time=slice('1-1-1991', '31-12-2020')).mean(dim='ensemble_member').mean(dim='time')
pr_carib_dry = pr_carib_dry.rolling(time=30*5,center=False).mean(dim='time').dropna("time")
pr_carib_dry_mean = pr_carib_dry.mean(dim='ensemble_member') 
pr_carib_dry_quant = pr_carib_dry.quantile([0.05,0.95],dim='ensemble_member') 

#------------------------OUTPUT------------------------#

fig, ax = plt.subplots(2,2 ,figsize=(16,12))

# precipitation
ax[0,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,(pr_carib_wet_mean-pr_carib_wet_ref_mean)*(100/pr_carib_wet_ref_mean),'steelblue',label='ave')
# ax[0,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,(pr_carib_wet_quant[0,:]-pr_carib_wet_ref_mean)*(100/pr_carib_wet_ref_mean),'blue',label='P05')
# ax[0,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,(pr_carib_wet_quant[1,:]-pr_carib_wet_ref_mean)*(100/pr_carib_wet_ref_mean),'red',label='P95')
ax[0,1].fill_between(tas_glob_wet_mean-tas_glob_wet_ref_mean,(pr_carib_wet_quant[0,:]-pr_carib_wet_ref_mean)*(100/pr_carib_wet_ref_mean),(pr_carib_wet_quant[1,:]-pr_carib_wet_ref_mean)*(100/pr_carib_wet_ref_mean),color='lightskyblue',linewidth=0)
ax[0,1].grid(True)
ax[0,1].set_ylim((-20,10))
ax[0,1].set_xlabel(r'$\Delta$ Global temperature rise $[^{\circ} C]$',fontsize = 16)
ax[0,1].set_ylabel(r'$\Delta$ Precipitation response $[\%]$',fontsize = 16)
ax[0,1].set_title('Wet season (May-November)',loc = 'right',fontsize = 14) 

ax[0,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,(pr_carib_dry_mean-pr_carib_dry_ref_mean)*(100/pr_carib_dry_ref_mean),'steelblue',label='ave')
# ax[0,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,(pr_carib_dry_quant[0,:]-pr_carib_dry_ref_mean)*(100/pr_carib_dry_ref_mean),'blue',label='P05')
# ax[0,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,(pr_carib_dry_quant[1,:]-pr_carib_dry_ref_mean)*(100/pr_carib_dry_ref_mean),'red',label='P95')
ax[0,0].fill_between(tas_glob_dry_mean-tas_glob_dry_ref_mean,(pr_carib_dry_quant[0,:]-pr_carib_dry_ref_mean)*(100/pr_carib_dry_ref_mean),(pr_carib_dry_quant[1,:]-pr_carib_dry_ref_mean)*(100/pr_carib_dry_ref_mean),color='lightskyblue',linewidth=0)
ax[0,0].set_ylim((-20,10))
ax[0,0].grid(True)
ax[0,0].set_xlabel(r'$\Delta$ Global temperature rise $[^{\circ} C]$',fontsize = 16)
ax[0,0].set_ylabel(r'$\Delta$ Precipitation response $[\%]$',fontsize = 16)
ax[0,0].set_title('Dry season (December-April)',loc = 'right',fontsize = 14)

percentiles =  mpatches.Patch( color='lightskyblue', label=r'$90\%$ EC-Earth3 ')
mean = mlines.Line2D([], [], color='steelblue', label='EC-Earth3 mean')
ax[0,1].legend(handles=[percentiles,mean],fontsize=14)

# temperature
ax[1,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,tas_carib_wet_mean-tas_carib_wet_ref_mean,'saddlebrown',label='ave')
# ax[1,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,tas_carib_wet_quant[0,:]-tas_carib_wet_ref_mean,'blue',label='P05')
# ax[1,1].plot(tas_glob_wet_mean-tas_glob_wet_ref_mean,tas_carib_wet_quant[1,:]-tas_carib_wet_ref_mean,'red',label='P95')
ax[1,1].fill_between(tas_glob_wet_mean-tas_glob_wet_ref_mean,tas_carib_wet_quant[0,:]-tas_carib_wet_ref_mean,tas_carib_wet_quant[1,:]-tas_carib_wet_ref_mean,color='sandybrown',linewidth=0)
ax[1,1].grid(True)
ax[1,1].set_ylim((-1,5))
ax[1,1].set_xlabel(r'$\Delta$ Global temperature rise $[^{\circ} C]$',fontsize = 16)
ax[1,1].set_ylabel(r'$\Delta$ Temperature response $^{\circ} C$',fontsize = 16)
ax[1,1].set_title('Wet season (May-November)',loc = 'right',fontsize = 14) 

ax[1,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,tas_carib_dry_mean-tas_carib_dry_ref_mean,'saddlebrown',label='ave')
# ax[1,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,tas_carib_dry_quant[0,:]-tas_carib_dry_ref_mean,'blue',label='P05')
# ax[1,0].plot(tas_glob_dry_mean-tas_glob_dry_ref_mean,tas_carib_dry_quant[1,:]-tas_carib_dry_ref_mean,'red',label='P95')
ax[1,0].fill_between(tas_glob_dry_mean-tas_glob_dry_ref_mean,tas_carib_dry_quant[0,:]-tas_carib_dry_ref_mean,tas_carib_dry_quant[1,:]-tas_carib_dry_ref_mean,color='sandybrown',linewidth=0)
ax[1,0].grid(True)
ax[1,0].set_ylim((-1,5))
ax[1,0].set_xlabel(r'$\Delta$ Global temperature rise $[^{\circ} C]$',fontsize = 16)
ax[1,0].set_ylabel(r'$\Delta$ Temperature response $^{\circ} C$',fontsize = 16)
ax[1,0].set_title('Dry season (December-April)',loc = 'right',fontsize = 14) 

percentiles =  mpatches.Patch( color='sandybrown', label=r'$90\%$ EC-Earth3 ')
mean = mlines.Line2D([], [], color='saddlebrown', label='EC-Earth3 mean')
ax[1,1].legend(handles=[percentiles,mean],fontsize=14)

ax[1,1].tick_params(axis='both', which='major', labelsize=14) # set all the x and y ticks to the same size
ax[1,0].tick_params(axis='both', which='major', labelsize=14)
ax[0,1].tick_params(axis='both', which='major', labelsize=14)
ax[0,0].tick_params(axis='both', which='major', labelsize=14)


fig.suptitle('Change EC-Earth3 ensembles statistics (1999-2110)',fontsize=20)
fig.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/definitivos/BES_fig1_Lenderink2014.png')

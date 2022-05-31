
# Python script to plot a map of the Atlantic ocean basin + Caribbean region + BES islands. (Figure 2 report)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.lines as mlines

plt.figure()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([0, 250, 0, 40])



ax.stock_img() # Put a background image on for nice sea rendering.

ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='gray'))
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='face', facecolor='black'))
ax.plot([-50, -50],[5,20],color = 'black',linewidth='.8')
ax.plot([-90,-90],[5,20],color = 'black',linewidth='.8')
ax.plot([-50,-90],[5,5],color = 'black',linewidth='.8')
ax.plot([-50,-90],[20,20],color = 'black',linewidth='.8')

ax.plot(-68.23,12.2, marker = 'o', color = 'black',markersize = 4 ,label = 'Bonaire')
ax.plot(-62.96,17.49, marker = 'o', color = 'red',markersize = 4, label = 'St.Eustatius & Saba' )
breakpoint()
Bonaire = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=4, label='Bonaire')
Saba = mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=4, label='St.Eustatius & Saba')
ax.legend(handles= [Bonaire, Saba], loc = 'upper left',fontsize = 10)

gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude= 0), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.ylines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 9, 'color': 'black'}
gl.ylabel_style = {'size': 9, 'color': 'black'}

plt.savefig('/home/bbrotonsbl/shared_data/volume_2/bbrotonsbl/resampling/Method-2014/output_plots/report/Caribbean_map.png', dpi = 1200)

import netCDF4
import numpy as np
from numpy import arange, meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
ncfile = netCDF4.Dataset('/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/EHF_heatwaves____yearly_summer.nc')
lon = ncfile.variables['lon']
lat = ncfile.variables['lat']
hwf = ncfile.variables['HWF_EHF'][:]
early = np.mean(hwf[:30,...], axis=0)
ave = np.mean(hwf, axis=0)
late = np.mean(hwf[-30:,...], axis=0)
parallels = arange(-40., -9., 10.)
meridians = arange(120., 160., 10.)
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4.5))
count = 0
colmap = 'YlOrRd'
for ax in axes.flat:
    map_axes = Basemap(ax=ax,projection='cyl',llcrnrlat=-44,urcrnrlat=-10,llcrnrlon=112, urcrnrlon=156)
    x, y = map_axes(*meshgrid(lon, lat))
    if count==0: 
        color = map_axes.pcolormesh(x,y,early,cmap=colmap)
    elif count==1:
        color = map_axes.pcolormesh(x,y,late,cmap=colmap)
    elif count==2: 
        color = map_axes.pcolormesh(x,y,ave,cmap=colmap)
    color.set_clim(vmin=0,vmax=30)
    map_axes.drawcoastlines()
    if count==0:
        map_axes.drawparallels(parallels,labels=[True,False,False,True],linewidth=0.0)
    else:
        map_axes.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.0)
    map_axes.drawmeridians(meridians,labels=[True,False,False,True],linewidth=0.0)
    count += 1
cbar_ax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat], orientation='horizontal')
cb = plt.colorbar(color, cax=cbar_ax, orientation='horizontal')
#plt.figtext(0.05,0.55,'a)')
plt.title('Mean HWF (Days)')
axes[0].set_title('a) 1911-1941')
axes[1].set_title('b) 1984-2013')
axes[2].set_title('c) 1911-2013')
#plt.show()
plt.savefig('HWF_basic_plots.eps',format='eps')

hwf = ncfile.variables['HWA_EHF'][:]
early = np.mean(hwf[:30,...], axis=0)
ave = np.mean(hwf, axis=0)
late = np.mean(hwf[-30:,...], axis=0)
parallels = arange(-40., -9., 10.)
meridians = arange(120., 160., 10.)
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4.5))
count = 0
colmap = 'YlOrRd'
missing = np.logical_not(early.mask==late.mask)
missing = missing/2.
missing[0,0] = 1
for ax in axes.flat:
    map_axes = Basemap(ax=ax,projection='cyl',llcrnrlat=-44,urcrnrlat=-10,llcrnrlon=112, urcrnrlon=156)
    x, y = map_axes(*meshgrid(lon, lat))
    if count==0:
        misscol = map_axes.pcolormesh(x,y,missing,cmap='Greys')
        color = map_axes.pcolormesh(x,y,early,cmap=colmap)
    elif count==1:
        color = map_axes.pcolormesh(x,y,late,cmap=colmap)
    elif count==2: 
        color = map_axes.pcolormesh(x,y,ave,cmap=colmap)
    color.set_clim(vmin=0,vmax=50)
    map_axes.drawcoastlines()
    if count==0:
        map_axes.drawparallels(parallels,labels=[True,False,False,True],linewidth=0.0)
    else:
        map_axes.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.0)
    map_axes.drawmeridians(meridians,labels=[True,False,False,True],linewidth=0.0)
    count += 1
cbar_ax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat], orientation='horizontal')
cb = plt.colorbar(color, cax=cbar_ax, orientation='horizontal')
#plt.figtext(0.05,0.55,'b)')
plt.title('Mean HWA ($^\circ C^2$)')
axes[0].set_title('d) 1911-1941')
axes[1].set_title('e) 1984-2013')
axes[2].set_title('f) 1911-2013')
#plt.show()
plt.savefig('HWA_basic_plots.eps',format='eps')

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:35:03 2017

@author: Tammas Loughran
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as bm

# Load the forcing data
ninofile = '/srv/ccrc/data35/z5032520/ancilforge/modal_anomalies.nc'
ninonc = nc.Dataset(ninofile,'r')
nino = ninonc.variables['elnino_sst'][:]
ninopic = nino[11:14,...].mean(axis=0)
lons2 = ninonc.variables['lon'][:]
lats2 = ninonc.variables['lat'][:]
nina = ninonc.variables['lanina_sst'][:]
ninapic = nina[11:14,...].mean(axis=0)

# Define the Figure Layout
fig, axes = plt.subplots(ncols=2,nrows=2,figsize=(9,6))

# Plot the El nino map
m = Basemap(ax=axes[0,0],projection='robin', lon_0=180.)
m.drawcoastlines()
parallels = np.arange(-90., 120., 30.)
meridians = np.arange(0., 360., 30.)
par = m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
for pr in par:
    try:
        par[pr][1][0].set_fontsize(9)
    except: pass
meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
for mr in meridians:
    try: 
        meridians[mr][1][0].set_rotation(50)
        meridians[mr][1][0].set_fontsize(9)
    except: pass
ninopic, lons = bm.shiftgrid(1, ninopic, lons2)
xx,yy = np.meshgrid(lons, lats2)
x,y = m(xx,yy)
fcont = m.contourf(x,y,ninopic,cmap='seismic',levels=np.arange(-2,2.1,0.2))
m.fillcontinents(color='gray',lake_color='gray')
cbar = m.colorbar(fcont, location='right', pad=0.3)
plt.title('a)', loc='left')

# Plot the lanina map
m = Basemap(ax=axes[0,1], projection='robin', lon_0=180.)
m.drawcoastlines()
par = m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
for pr in par:
    try:
        par[pr][1][0].set_fontsize(9)
    except: pass
meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
for mr in meridians:
    try: 
        meridians[mr][1][0].set_rotation(50)
        meridians[mr][1][0].set_fontsize(9)
    except: pass
ninapic, _ = bm.shiftgrid(1, ninapic, lons2)
xx,yy = np.meshgrid(lons, lats2)
x,y = m(xx,yy)
fcont = m.contourf(x,y,ninapic,cmap='seismic',levels=np.arange(-2,2.1,0.2))
m.fillcontinents(color='gray',lake_color='gray')
cbar = m.colorbar(fcont, location='right', pad=0.3)
plt.title('b)', loc='left')

# Plot the time series
i, j = 189, 96
nno_series = nino[:,j,i]
nna_series = nina[:,j,i]

plt.sca(axes[1,0])
plt.plot(nno_series,'k-',label='El Nino')
plt.scatter(np.arange(0,24), nno_series, marker='x',color='k')
plt.plot(nna_series,'k',label='La Nina', linestyle='dotted')
plt.scatter(np.arange(0,24), nna_series, marker='x',color='k')
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels.extend(labels)
plt.xticks(np.arange(0,24), labels, rotation=45, fontsize=9)
plt.axhline(y=0,color='k')
plt.xlim(0, 24)
plt.grid(True)
plt.xlabel('Month')
plt.ylabel('SST $^{\circ}$C')
plt.legend(loc='upper right')
plt.title('c)', loc='left')


# Plot the box regions.
def plotbox_onmap(lat1,lat2,lon1,lon2,bm):
    x,y = m([lon1,lon2],[lat1,lat2])
    m.plot([x[0],x[0]],[y[0],y[1]],color='k')
    m.plot([x[0],x[1]],[y[1],y[1]],color='k')
    m.plot([x[1],x[1]],[y[1],y[0]],color='k')
    m.plot([x[1],x[0]],[y[0],y[0]],color='k')
    m.drawcoastlines()


m = Basemap(ax=axes[1,1],projection='mill',
            llcrnrlon=112.5,llcrnrlat=-43.75,
            urcrnrlon=155.625,urcrnrlat=-10)
#regions = [(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)]
regions = [(129.,-12),(139.,-18),(145.5,-24),(141,-31)]
m.drawmeridians([120,130,140,150],labels=[True,False,False,True],fontsize=9)
m.drawparallels([10,20,30,40],labels=[True,False,False,True],fontsize=9)
for i in xrange(len(regions)):
    plotbox_onmap(regions[i][1],regions[i][1]-7.5,regions[i][0],regions[i][0]+7.5,m)
plt.sca(axes[1,1])
plt.title('d)', loc='left')
plt.subplots_adjust(hspace=0.4,wspace=0.3)
plt.savefig('ENSO_forcing.eps',format='eps')

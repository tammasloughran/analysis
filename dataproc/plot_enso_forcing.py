# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:35:03 2017

@author: Tammas Loughran
"""

import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as bm
from scipy import signal
import sys

ninofile = '/srv/ccrc/data35/z5032520/ancilforge/modal_anomalies.nc'
ninonc = nc.Dataset(ninofile,'r')
nino = ninonc.variables['elnino_sst'][:]
ninopic = nino[11:14,...].mean(axis=0)
lons2 = ninonc.variables['lon'][:]
lats2 = ninonc.variables['lat'][:]
nina = ninonc.variables['lanina_sst'][:]
ninapic = nina[11:14,...].mean(axis=0)

fig, axes = plt.subplots(nrows=3,figsize=(7,10))

m = Basemap(ax=axes[0],projection='robin', lon_0=180.)
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
#m.contour(x,y,modokipic,colors='k',levels=[0])
m.fillcontinents(color='gray',lake_color='gray')
#cbar = plt.colorbar(fcont, orientation='horizontal')
cbar = m.colorbar(fcont, location='right', pad=0.3)
#cbar.ax.get_xaxis().set_ticks([])
#for j, lab in enumerate([str(i) for i in np.arange(-1.5,1.7,0.2)]):
#    cbar.ax.text((1/15.)*j, -0.6, lab, ha='center', va='center')
plt.title('a)', loc='left')

m = Basemap(ax=axes[1], projection='robin', lon_0=180.)
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

i, j = 189, 96
nno_series = nino[:,j,i]
nna_series = nina[:,j,i]

plt.sca(axes[2])
plt.plot(nno_series,'k-',label='El Nino')
plt.scatter(np.arange(0,24), nno_series, marker='x',color='k')
plt.plot(nna_series,'k',label='La Nina', linestyle='dotted')
plt.scatter(np.arange(0,24), nna_series, marker='x',color='k')
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels.extend(labels)
plt.xticks(np.arange(0,24), labels, rotation=45)
plt.axhline(y=0,color='k')
plt.xlim(0, 24)
plt.grid(True)
plt.xlabel('Month')
plt.ylabel('SST $^{\circ}$C')
plt.legend(loc='upper right')
plt.title('c)', loc='left')
plt.savefig('ENSO_forcing.eps',format='eps')
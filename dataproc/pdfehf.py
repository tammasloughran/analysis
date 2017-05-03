# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:13:10 2017

@author: Tammas Loughran
"""

import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 

directory = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
model = nc.Dataset(directory+'vamrb/EHF_heatwaves_ACCESS1.3_vamrb_yearly_summer.nc')
model_tx90 = model.variables['t90pct'][:]

directory = '/srv/ccrc/data35/z5032520/ehfheatwaves/'
obs = nc.Dataset(directory+'ehf_remap.nc')
obs_tx90 = obs.variables['t90pct'][:] 
lons = obs.variables['lon'][:]
lats = obs.variables['lat'][:]
ausx = (lons>110)&(lons<155)
ausy = (lats>-45)&(lats<-10)
aus_lons = lons[ausx]
aus_lats = lats[ausy]


dates = pd.date_range(start='2001-01-01',end='2001-12-31',freq='D')
summer = (dates.month==11)|(dates.month==12)|(dates.month==1)|(dates.month==2)|(dates.month==3)

diff = model_tx90 - obs_tx90

diff_aus = diff[:,:,ausx][:,ausy,:]

dsummer = diff_aus[summer,...].mean(axis=0)

mp = Basemap(projection='cyl',
        llcrnrlat=-45, urcrnrlat=-10,
        llcrnrlon=110, urcrnrlon=155,resolution='l')
x,y = np.meshgrid(aus_lons-0.5,aus_lats-0.5)
i,j = np.meshgrid(aus_lons,aus_lats)
ii,jj = mp(i,j)
xx,yy = mp(x,y)
colour = mp.pcolormesh(xx,yy,dsummer,vmin=-6,vmax=6,cmap='bwr')
cont = mp.contour(ii,jj,dsummer,levels=range(-6,7,1),colors='k',linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
cb = mp.colorbar(colour,ticks=range(-6,7,1))
cb.ax.set_xlabel('$^\circ$C')
mp.drawcoastlines()
plt.savefig('tx90pct_diff_obs_access.png',format='png')
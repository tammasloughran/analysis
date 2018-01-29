# -*- coding: utf-8 -*-
"""
Created on Wed January  17 14:45:42 2017

@author: Tammas Loughran
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats


# Load the heatwave file
hwfile = '/srv/ccrc/data48/z5032520/20crv2/EHF_heatwaves_20crv2_1901-2012_daily.nc'
hwnc = nc.Dataset(hwfile,'r')
hwevents = hwnc.variables['event'][:] 
hwlats = hwnc.variables['lat'][:]
hwlons = hwnc.variables['lon'][:]
hwdates = nc.num2date(hwnc.variables['time'][:],hwnc.variables['time'].units)
hwyears = np.array([d.year for d in hwdates])
hwmonths = np.array([d.month for d in hwdates])
hwnc.close()

# Load the land sea mask
obsdir = '/srv/ccrc/data48/z5032520/20crv2/'
lsmnc = nc.Dataset(obsdir+'20CRV2c_mask.nc', 'r')
lsm = lsmnc.variables['lsm'][:]
lsmnc.close()

# Define El Nino and La Nina years. year is year containing December
elninoyears = [1911,1913,1914,1918,1925,1930,1941,1951,1957,1965,1969,1972,1976,1982,1987,1997,2006]
laninayears = [1909,1915,1920,1933,1942,1949,1950,1955,1970,1973,1975,1984,1988,1998,1999,2007,2010]
#elninoyears = [1972,1976,1982,1987,1997,2006]
#laninayears = [1970,1973,1975,1984,1988,1998,1999,2007,2010]

startyear = 1901

# Load the lhtfl data
lhtfldir = '/srv/ccrc/data48/z5032520/20crv2/lhtfl'
ncfile = nc.MFDataset(lhtfldir+'/lhtfl.*.nc','r')
dates = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units)
years = np.array([i.year for i in dates])
lhtfl = ncfile.variables['lhtfl'][years>=startyear,...]
dates = dates[years>=startyear]
years = years[years>=startyear]
months = np.array([d.month for d in dates])
novmar = (months==11)|(months==12)|(months==1)|(months==2)|(months==3)
lhtflclim = lhtfl[novmar,...].mean(axis=0)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Load the shtfl data
shtfldir = '/srv/ccrc/data48/z5032520/20crv2/shtfl'
ncfile = nc.MFDataset(shtfldir+'/shtfl.*.nc','r')
dates = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units)
years = np.array([i.year for i in dates])
shtfl = ncfile.variables['shtfl'][years>=startyear,...]
dates = dates[years>=startyear]
years = years[years>=startyear]
months = np.array([d.month for d in dates])
novmar = (months==11)|(months==12)|(months==1)|(months==2)|(months==3)
shtflclim = shtfl[novmar,...].mean(axis=0)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Define index of nov 1st and mar 31
nov1 = 304
mar31 = 454

shape = (0,)+(ncfile.variables['shtfl'].shape[1:])
snino = np.ones(shape)
lnino = np.ones(shape)

for year in elninoyears:
    print(year)
    sflyear = shtfl[(years==year)|(years==year+1),...]
    if year%4==0: sflyear = np.delete(sflyear,59,axis=0)
    if (year+1)%4==0: sflyear = np.delete(sflyear,424,axis=0)
    sflyear = sflyear[nov1:mar31+1,...]
    lflyear = lhtfl[(years==year)|(years==year+1),...]
    if year%4==0: lflyear = np.delete(lflyear,59,axis=0)
    if (year+1)%4==0: lflyear = np.delete(lflyear,424,axis=0)
    lflyear = lflyear[nov1:mar31+1,...]
    snino = np.append(snino, sflyear.mean(axis=0)[None,...], axis=0)
    lnino = np.append(lnino, lflyear.mean(axis=0)[None,...], axis=0)

snina = np.ones(shape)
lnina = np.ones(shape)

print('La Nina')
for year in laninayears:
    print(year)
    sflyear = shtfl[(years==year)|(years==year+1),...]
    if year%4==0: sflyear = np.delete(sflyear,59,axis=0)
    if (year+1)%4==0: sflyear = np.delete(sflyear,424,axis=0)
    sflyear = sflyear[nov1:mar31+1,...]
    lflyear = lhtfl[(years==year)|(years==year+1),...]
    if year%4==0: lflyear = np.delete(lflyear,59,axis=0)
    if (year+1)%4==0: lflyear = np.delete(lflyear,424,axis=0)
    lflyear = lflyear[nov1:mar31+1,...]
    snina = np.append(snina, sflyear.mean(axis=0)[None,...], axis=0)
    lnina = np.append(lnina, lflyear.mean(axis=0)[None,...], axis=0)

sdiff = snino.mean(axis=0)-snina.mean(axis=0)
ldiff = lnino.mean(axis=0)-lnina.mean(axis=0)


def plot_Q(data,ax,levels,ndays=0,colormap='YlOrRd'):
    data = np.ma.array(data, mask=np.logical_not(lsm))
    m = Basemap(ax=ax,projection='mill',
                llcrnrlon=105.,llcrnrlat=-45.,
                urcrnrlon=160.,urcrnrlat=-5.)
    lns,lts = np.meshgrid(lons-1,lats+1)
    mshx,mshy = m(lns,lts)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    cont = m.pcolormesh(mshx,mshy,data,cmap=colormap,vmin=levels[0],vmax=levels[-1])
    cnt = m.contour(x,y,data,colors='k',linewidths=0.3,levels=levels)
    for c in cnt.collections:
        if c.get_linestyle() == [(None, None)]: continue
        else: c.set_dashes([(0, (2.0, 2.0))])
    m.drawcoastlines()
    cbar = m.colorbar(cont,location='right',ticks=levels)
    cbar.set_label('$Wm^{-2}$')
    return cont, m


f, axes = plt.subplots(nrows=2, ncols=2,figsize=(7,6.5))

_, m = plot_Q(shtflclim, axes[0][0], np.arange(0,161,20), colormap='YlOrRd')
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
axes[0][0].set_title('a) Q$_{H}$ Climatology', loc='left')

_, m = plot_Q(lhtflclim, axes[0][1], np.arange(0,135,15), colormap='YlGnBu')
axes[0][1].set_title('b) Q$_{E}$ Climatology', loc='left')

_, m = plot_Q(sdiff, axes[1][0], np.arange(-30,31,10), colormap='bwr')
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
axes[1][0].set_title('c) Q$_{H}$ El Nino - La Nina', loc='left')

_, m = plot_Q(ldiff, axes[1][1], np.arange(-20,21,5), colormap='BrBG')
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
axes[1][1].set_title('d) Q$_{E}$ El Nino - La Nina', loc='left')

plt.subplots_adjust(hspace=0.1,wspace=0.3)
#plt.show()
plt.savefig('Q_alldays.eps',format='eps')




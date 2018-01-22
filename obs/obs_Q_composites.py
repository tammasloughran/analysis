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
#elninoyears = [1911,1913,1914,1918,1925,1930,1941,1951,1957,1965,1969,1972,1976,1982,1987,1997,2006]
#laninayears = [1909,1915,1920,1933,1942,1949,1950,1955,1970,1973,1975,1984,1988,1998,1999,2007,2010]
elninoyears = [1972,1976,1982,1987,1997,2006]
laninayears = [1970,1973,1975,1984,1988,1998,1999,2007,2010]


# Load the lhtfl data
lhtfldir = '/srv/ccrc/data48/z5032520/20crv2/lhtfl'
ncfile = nc.MFDataset(lhtfldir+'/lhtfl.*.nc','r')
dates = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units)
years = np.array([i.year for i in dates])
lhtfl = ncfile.variables['lhtfl'][years>=1970,...]
dates = dates[years>=1970]
years = years[years>=1970]
months = np.array([d.month for d in dates])
novmar = (months==11)|(months==12)|(months==1)|(months==2)|(months==3)
clim = lhtfl[novmar,...].mean(axis=0)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Define index of nov 1st and mar 31
nov1 = 304
mar31 = 454

# Define finction to select region from hwdata
def select_region(lon, lat):
    x = (hwlons>=lon)&(hwlons<=lon+7.5)
    y = (hwlats<=lat)&(hwlats>=lat-7.5)
    region = hwdata[:,y,:][:,:,x]
    region = region.sum(axis=1).sum(axis=1)
    region = region>=9
    return region

# Note: region top left coordinates
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)

# Create el nino arrays
shape = (0,)+(ncfile.variables['lhtfl'].shape[1:])
seaus_o = np.ones(shape)
neaus_o = seaus_o.copy()
naus_o = seaus_o.copy()
eaus_o = seaus_o.copy()

print('Looping el nino')
#loop over el nino years
for year in elninoyears:
    print(year)
    # Fetch the Q data
    htfl = lhtfl[(years==year)|(years==year+1)]
    if year%4==0: # is a leap year
        htfl = np.delete(htfl,59,axis=0)
    if (year+1)%4==0:
        htfl = np.delete(htfl,424,axis=0)
    htfl = htfl[nov1:mar31+1,...]
    # Select the heatwave data for this year
    hwdata = hwevents[(hwyears==year)|(hwyears==year+1),...]
    hwdata = hwdata[nov1:mar31+1,...]
    # Select the region of interest
    seaus = select_region(141.,-31.)
    neaus = select_region(139.,-18.)
    naus = select_region(129.,-12)
    eaus = select_region(145.5,-24)
    # Add to array
    seaus_o = np.append(seaus_o, htfl[seaus,...],axis=0)
    neaus_o = np.append(neaus_o, htfl[neaus,...],axis=0)
    naus_o = np.append(naus_o, htfl[naus,...],axis=0)
    eaus_o = np.append(eaus_o, htfl[eaus,...],axis=0)

# Create la nina arrays
shape = (0,)+(ncfile.variables['lhtfl'].shape[1:])
seaus_a = np.ones(shape)
neaus_a = seaus_a.copy()
naus_a = seaus_a.copy()
eaus_a = seaus_a.copy()

print('Looping la nina')
# Loop over all la nina years
for year in laninayears:
    print(year)
    # Load the mslp data
    htfl = lhtfl[(years==year)|(years==year+1)]
    if year%4==0: # is a leap year
        htfl = np.delete(htfl,59,axis=0)
    if (year+1)%4==0:
        htfl = np.delete(htfl,424,axis=0)
    htfl = htfl[nov1:mar31+1,...]
    # Select the heatwave data for this year
    hwdata = hwevents[(hwyears==year)|(hwyears==year+1),...]
    hwdata = hwdata[nov1:mar31+1,...]
    # Select the region of interest
    seaus = select_region(141.,-31.)
    neaus = select_region(139.,-18.)
    naus = select_region(129.,-12)
    eaus = select_region(145.5,-24)
    # Add to array
    seaus_a = np.append(seaus_a, htfl[seaus,...],axis=0)
    neaus_a = np.append(neaus_a, htfl[neaus,...],axis=0)
    naus_a = np.append(naus_a, htfl[naus,...],axis=0)
    eaus_a = np.append(eaus_a, htfl[eaus,...],axis=0)


def plot_Q(data,ax,ndays=0):
    data = np.ma.array(data, mask=np.logical_not(lsm))
    m = Basemap(ax=ax,projection='mill',
                llcrnrlon=105.,llcrnrlat=-45.,
                urcrnrlon=160.,urcrnrlat=-5.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    levels = np.arange(-50,51,10)
    cont = m.pcolormesh(x,y,data,cmap='PuOr',vmin=-50,vmax=50)
    cnt = m.contour(x,y,data,colors='k',linewidths=0.3,levels=levels)
    for c in cnt.collections:
        if c.get_linestyle() == [(None, None)]: continue
        else: c.set_dashes([(0, (2.0, 2.0))])
    m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
    m.drawcoastlines()
    return cont, m

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(6,7.75))
# El nino
# Express as anomaly ignificance and plot 
data = seaus_o.mean(axis=0) - clim
t, sig = stats.ttest_ind(seaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[3][0])
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
ndays = seaus_o.shape[0]
axes[3][0].set_title('g) SE n='+str(ndays), loc='left')

data = eaus_o.mean(axis=0) - clim
t, sig = stats.ttest_ind(eaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[2][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = eaus_o.shape[0]
axes[2][0].set_title('e) E n='+str(ndays), loc='left')

data = neaus_o.mean(axis=0) - clim
t, sig = stats.ttest_ind(neaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[1][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = neaus_o.shape[0]
axes[1][0].set_title('c) NE n='+str(ndays), loc='left')

data = naus_o.mean(axis=0) - clim
t, sig = stats.ttest_ind(naus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[0][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = naus_o.shape[0]
axes[0][0].set_title('a) N n='+str(ndays), loc='left')


# La Nina
data = seaus_a.mean(axis=0) - clim
t, sig = stats.ttest_ind(seaus_a,lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[3][1])
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = seaus_a.shape[0]
axes[3][1].set_title('h) SE n='+str(ndays), loc='left')

data = eaus_a.mean(axis=0) - clim
t, sig = stats.ttest_ind(eaus_a,lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[2][1])
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = eaus_a.shape[0]
axes[2][1].set_title('f) E n='+str(ndays), loc='left')

data = neaus_a.mean(axis=0) - clim
t, sig = stats.ttest_ind(neaus_a, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[1][1])
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = neaus_a.shape[0]
axes[1][1].set_title('d) NE n='+str(ndays), loc='left')

data = naus_a.mean(axis=0) - clim
t, sig = stats.ttest_ind(naus_a, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
contours, m = plot_Q(data,axes[0][1])
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = naus_a.shape[0]
axes[0][1].set_title('b) N n='+str(ndays), loc='left')


cax = f.add_axes([0.1,0.07,0.8,0.02])
cbar = plt.colorbar(contours,cax=cax,orientation='horizontal',ticks=np.arange(-50,51,10))
cbar.set_label('$Wm^{-2}$')
f.suptitle('El Nino            La Nina', fontsize=20)
plt.savefig('lhtfl1970_hwdays_composites.eps',format='eps')
#plt.show()
print 'done'

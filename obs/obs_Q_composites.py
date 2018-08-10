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
import pdb

## USER
vname = 'lhtfl'

# Load the heatwave file
print('Loading data')
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


# Load the lhtfl data
lhtfldir = '/srv/ccrc/data48/z5032520/20crv2/'+vname
ncfile = nc.MFDataset(lhtfldir+'/'+vname+'.*.nc','r')
dates = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units)
years = np.array([i.year for i in dates])
lhtfl = ncfile.variables[vname][...]
#lhtfl = ncfile.variables['lhtfl'][years>=1970,...]
#dates = dates[years>=1970]
#years = years[years>=1970]
months = np.array([d.month for d in dates])
novmar = (months==11)|(months==12)|(months==1)|(months==2)|(months==3)
# The climatology should be subtracted from the data for each month.
clim1 = np.ones((5,)+lhtfl.shape[1:])*np.nan
for i,m in enumerate([11,12,1,2,3]):
    clim1[i,...] = lhtfl[months==m,...].mean(axis=0)
clim = np.ones((151,)+lhtfl.shape[1:])*np.nan
clim[0:30,...] = clim1[0]
clim[30:61,...] = clim1[1]
clim[61:92,...] = clim1[2]
clim[92:120,...] = clim1[3]
clim[120:151,...] = clim1[4] # This constructs a climatology to be subtracted from a nov-mar chunck of data.
#clim = lhtfl[novmar,...].mean(axis=0) # This was the old climatology. I think it's shit.
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]



# Define index of nov 1st and mar 31
nov1 = 304
mar31 = 454

# E(61,79), NE(58,78), SE(65,78), N(55,74)
east = (61,79)
neast = (58,78)
seast = (65,78)
print('east',lats[east[0]],'N',lons[east[1]],'E')
print('northeast',lats[neast[0]],'N',lons[neast[1]],'E')
print('southeast',lats[seast[0]],'N',lons[seast[1]],'E')
yy,xx = 58,78

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
shape = (0,)+(ncfile.variables[vname].shape[1:])
seaus_o = np.ones(shape)*np.nan
neaus_o = seaus_o.copy()
naus_o = seaus_o.copy()
eaus_o = seaus_o.copy()
obstso = np.ones((0,))*np.nan

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
    htfl = htfl[nov1:mar31+1,...] - clim
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
    if not seaus.any():
        obstso = np.append(obstso, np.nan)
    else:
        obstso = np.append(obstso, htfl[seaus,yy,xx].mean(axis=0))


# Create la nina arrays
shape = (0,)+(ncfile.variables[vname].shape[1:])
seaus_a = np.ones(shape)*np.nan
neaus_a = seaus_a.copy()
naus_a = seaus_a.copy()
eaus_a = seaus_a.copy()
obstsa = np.ones((0,))

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
    htfl = htfl[nov1:mar31+1,...] - clim
    # Select the heatwave data for this year
    hwdata = hwevents[(hwyears==year)|(hwyears==year+1),...]
    hwdata = hwdata[nov1:mar31+1,...]
    # Select the region of interest
    seaus = select_region(141.,-31.)
    neaus = select_region(139.,-18.)
    naus = select_region(129.,-12)
    eaus = select_region(145.5,-24)
    #print(np.where(naus>0))
    # Add to array
    seaus_a = np.append(seaus_a, htfl[seaus,...],axis=0)
    neaus_a = np.append(neaus_a, htfl[neaus,...],axis=0)
    naus_a = np.append(naus_a, htfl[naus,...],axis=0)
    eaus_a = np.append(eaus_a, htfl[eaus,...],axis=0)
    if not seaus.any():
        obstsa = np.append(obstsa, np.nan)
    else:
        obstsa = np.append(obstsa, htfl[seaus,yy,xx].mean(axis=0))

def plot_Q(data,ax,ndays=0):
    if vname=='lhtfl':
        colors = 'PuOr'
    elif vname=='shtfl':
        colors = 'bwr'
    data = np.ma.array(data, mask=np.logical_not(lsm))
    m = Basemap(ax=ax,projection='mill',
                llcrnrlon=105.,llcrnrlat=-45.,
                urcrnrlon=160.,urcrnrlat=-5.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    levels = np.arange(-50,51,10)
    cont = m.pcolormesh(x,y,data,cmap=colors,vmin=-50,vmax=50)
    cnt = m.contour(x,y,data,colors='k',linewidths=0.3,levels=levels)
    for c in cnt.collections:
        if c.get_linestyle() == [(None, None)]: continue
        else: c.set_dashes([(0, (2.0, 2.0))])
    m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
    m.drawcoastlines()
    return cont, m

# This converts november to march months to anomalies
for i,m in enumerate([11,12,1,2,3]):
    lhtfl[months==m] = lhtfl[months==m] - clim1[i]

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(6,7.75))
# El nino
# Express as anomaly significance and plot 
data = seaus_o.mean(axis=0)# - clim
t, sig = stats.ttest_ind(seaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[3][0])
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
ndays = seaus_o.shape[0]
axes[3][0].set_title('g) SE n='+str(ndays), loc='left')

data = eaus_o.mean(axis=0)# - clim
t, sig = stats.ttest_ind(eaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[2][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = eaus_o.shape[0]
axes[2][0].set_title('e) E n='+str(ndays), loc='left')

data = neaus_o.mean(axis=0)# - clim
t, sig = stats.ttest_ind(neaus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[1][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = neaus_o.shape[0]
axes[1][0].set_title('c) NE n='+str(ndays), loc='left')

data = naus_o.mean(axis=0)# - clim
t, sig = stats.ttest_ind(naus_o, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[0][0])
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = naus_o.shape[0]
axes[0][0].set_title('a) N n='+str(ndays), loc='left')


# La Nina
data = seaus_a.mean(axis=0)# - clim
t, sig = stats.ttest_ind(seaus_a,lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[3][1])
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = seaus_a.shape[0]
axes[3][1].set_title('h) SE n='+str(ndays), loc='left')

data = eaus_a.mean(axis=0)# - clim
t, sig = stats.ttest_ind(eaus_a, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[2][1])
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = eaus_a.shape[0]
axes[2][1].set_title('f) E n='+str(ndays), loc='left')

data = neaus_a.mean(axis=0)# - clim
t, sig = stats.ttest_ind(neaus_a, lhtfl[novmar,...], axis=0, equal_var=False)
sig = np.ma.array(sig, mask=np.logical_not(lsm))
_, m = plot_Q(data,axes[1][1])
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = neaus_a.shape[0]
axes[1][1].set_title('d) NE n='+str(ndays), loc='left')

data = naus_a.mean(axis=0)# - clim
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
plt.savefig(vname+'_hwdays_composites.eps',format='eps')
#plt.show()

#Example gridpoint
plt.figure()
plt.axhline(np.nanmean(obstsa),color='b')
plt.axhline(np.nanmean(obstso),color='r')
plt.scatter(laninayears, obstsa)
plt.scatter(elninoyears, obstso)
plt.ylabel('$W/m^{2}$')
plt.xlabel('year')
plt.legend(['La Nina','El Nino'])
plt.savefig('gridpoint_scatter_'+vname+'_NE.eps',format='eps')

print 'done'

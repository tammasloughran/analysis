# -*- coding: utf-8 -*-
"""
Created on Wed January 17 14:45:42 2017

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
lats = hwnc.variables['lat'][:]
lons = hwnc.variables['lon'][:]
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


# Load the MSLP data
monthlydir = '/srv/ccrc/data35/z5032520/20CRv2/mslp'
monmslpnc = nc.MFDataset(monthlydir+'/monthly_prmsl.*.nc','r')
mmslp = monmslpnc.variables['prmsl'][:]
mlats = monmslpnc.variables['lat'][:]
mlons = monmslpnc.variables['lon'][:]
dates = nc.num2date(monmslpnc.variables['time'][:],monmslpnc.variables['time'].units)
years = np.array([i.year for i in dates])
#months = np.array(d.month for d in dates])
startyear = 1901
mmslp = mmslp[years>=startyear,...]
#clim = mmslp
elclim = np.ones((0,)+(mmslp.shape[1:]))
for year in elninoyears:
    index = ((year-startyear)*12)+10
    elclim = np.append(elclim, mmslp[index:index+5,...], axis=0)
elpclim = np.ones((5,)+elclim.shape[1:])*np.nan
for i in xrange(5):
    elpclim[i,...] = elclim[i::5].mean(axis=0)
elmslp_clim = elclim.mean(axis=0)

laclim = np.ones((0,)+(mmslp.shape[1:]))
for year in laninayears:
    index = ((year-startyear)*12)+10
    laclim = np.append(laclim, mmslp[index:index+5,...], axis=0)
lapclim = np.ones((5,)+laclim.shape[1:])*np.nan
for i in xrange(5):
    lapclim[i,...] = laclim[i::5].mean(axis=0)
lamslp_clim = laclim.mean(axis=0)

mslp_clim = mmslp.mean(axis=0)

# Load the ncfile for daily mslp data
files = [obsdir+'prmsl/prmsl.'+str(yr)+'.nc' for yr in range(startyear,2013)]
ncfile = nc.MFDataset(files)
dates = nc.num2date(ncfile.variables['time'][:],ncfile.variables['time'].units)
years = np.array([d.year for d in dates])
months = np.array([d.month for d in dates])

# Define index of nov 1st and mar 31
nov1 = 304
mar31 = 454

# Define finction to select region from hwdata
def select_region(lon, lat):
    x = (lons>=lon)&(lons<=lon+7.5)
    y = (lats<=lat)&(lats>=lat-7.5)
    region = hwdata[:,y,:][:,:,x]
    region = region.sum(axis=1).sum(axis=1)
    region = region>=9
    return region

# Note: region top left coordinates
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)

# Create el nino arrays
shape = (0,)+(ncfile.variables['prmsl'].shape[1:])
seaus_mslp_o = np.ones(shape)
neaus_mslp_o = seaus_mslp_o.copy()
naus_mslp_o = seaus_mslp_o.copy()
eaus_mslp_o = seaus_mslp_o.copy()

def anomalise(data,pclim):
    data[0:30,...] = data[0:30,...] - pclim[0,...]
    data[30:61,...] = data[30:61,...] - pclim[1,...]
    data[61:92,...] = data[61:92,...] - pclim[2,...]
    data[92:120,...] = data[92:120,...] - pclim[3,...]
    data[120:151,...] = data[120:151,...] - pclim[4,...]
    return data

print('Looping el nino')
#loop over el nino years
for year in elninoyears:
    print(year)
    # Load the mslp data
    mslpdata = ncfile.variables['prmsl'][(years==year)|(years==year+1)]
    if year%4==0: # is a leap year
        mslpdata = np.delete(mslpdata,59,axis=0)
    if (year+1)%4==0:
        mslpdata = np.delete(mslpdata,424,axis=0)
    mslpdata = mslpdata[nov1:mar31+1,...]
    mslpdata = anomalise(mslpdata,elpclim)
    # Select the heatwave data for this year
    hwdata = hwevents[(hwyears==year)|(hwyears==year+1),...]
    hwdata = hwdata[nov1:mar31+1,...]
    # Select the region of interest
    seaus = select_region(141.,-31.)
    neaus = select_region(139.,-18.)
    naus = select_region(129.,-12)
    eaus = select_region(145.5,-24)
    # Add to array
    seaus_mslp_o = np.append(seaus_mslp_o, mslpdata[seaus,...],axis=0)
    neaus_mslp_o = np.append(neaus_mslp_o, mslpdata[neaus,...],axis=0)
    naus_mslp_o = np.append(naus_mslp_o, mslpdata[naus,...],axis=0)
    eaus_mslp_o = np.append(eaus_mslp_o, mslpdata[eaus,...],axis=0)

# Create la nina arrays
shape = (0,)+(ncfile.variables['prmsl'].shape[1:])
seaus_mslp_a = np.ones(shape)
neaus_mslp_a = seaus_mslp_a.copy()
naus_mslp_a = seaus_mslp_a.copy()
eaus_mslp_a = seaus_mslp_a.copy()

print('Looping la nina')
# Loop over all la nina years
for year in laninayears:
    print(year)
    # Load the mslp data
    mslpdata = ncfile.variables['prmsl'][(years==year)|(years==year+1)]
    if year%4==0: # is a leap year
        mslpdata = np.delete(mslpdata,59,axis=0)
    if (year+1)%4==0:
        mslpdata = np.delete(mslpdata,424,axis=0)
    mslpdata = mslpdata[nov1:mar31+1,...]
    mslpdata = anomalise(mslpdata,lapclim)
    # Select the heatwave data for this year
    hwdata = hwevents[(hwyears==year)|(hwyears==year+1),...]
    hwdata = hwdata[nov1:mar31+1,...]
    # Select the region of interest
    seaus = select_region(141.,-31.)
    neaus = select_region(139.,-18.)
    naus = select_region(129.,-12)
    eaus = select_region(145.5,-24)
    # Add to array
    seaus_mslp_a = np.append(seaus_mslp_a, mslpdata[seaus,...],axis=0)
    neaus_mslp_a = np.append(neaus_mslp_a, mslpdata[neaus,...],axis=0)
    naus_mslp_a = np.append(naus_mslp_a, mslpdata[naus,...],axis=0)
    eaus_mslp_a = np.append(eaus_mslp_a, mslpdata[eaus,...],axis=0)


def plot_mslp(data,ax,ndays=0):
    m = Basemap(ax=ax,projection='mill',
                llcrnrlon=80.,llcrnrlat=-60.,
                urcrnrlon=220.,urcrnrlat=-5.)
    lns,lts = np.meshgrid(mlons,mlats)
    x,y = m(lns,lts)
    levels = np.arange(-700,800,100)
    cont = m.contourf(x,y,data,cmap='bwr', levels=levels, extend='both')
    m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=1)
    m.drawcoastlines()
    return cont, m

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(6,7.75))
# El nino
# Express as anomaly ignificance and plot 
data = seaus_mslp_o.mean(axis=0)
t, sig = stats.ttest_ind(seaus_mslp_o+elmslp_clim, elclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[3][0])
m.drawmeridians([100,140,180,220],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([0,-20,-40,-60],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
ndays = seaus_mslp_o.shape[0]
axes[3][0].set_title('g) SE n='+str(ndays), loc='left')

data = eaus_mslp_o.mean(axis=0)
t, sig = stats.ttest_ind(eaus_mslp_o+elmslp_clim, elclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[2][0])
m.drawparallels([0,-20,-40,-60],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([100,140,180,220],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = eaus_mslp_o.shape[0]
axes[2][0].set_title('e) E n='+str(ndays), loc='left')

data = neaus_mslp_o.mean(axis=0)
t, sig = stats.ttest_ind(neaus_mslp_o+elmslp_clim, elclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[1][0])
m.drawparallels([0,-20,-40,-60],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([100,140,180,220],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = neaus_mslp_o.shape[0]
axes[1][0].set_title('c) NE n='+str(ndays), loc='left')

data = naus_mslp_o.mean(axis=0)
t, sig = stats.ttest_ind(naus_mslp_o+elmslp_clim, elclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[0][0])
m.drawparallels([0,-20,-40,-60],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([100,140,180,220],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = naus_mslp_o.shape[0]
axes[0][0].set_title('a) N n='+str(ndays), loc='left')


# La Nina
data = seaus_mslp_a.mean(axis=0)
t, sig = stats.ttest_ind(seaus_mslp_a+lamslp_clim, laclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[3][1])
m.drawmeridians([100,140,180,220],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawparallels([0,-20,-40,-60],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ndays = seaus_mslp_a.shape[0]
axes[3][1].set_title('h) SE n='+str(ndays), loc='left')

data = eaus_mslp_a.mean(axis=0)
t, sig = stats.ttest_ind(eaus_mslp_a+lamslp_clim, laclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[2][1])
ndays = eaus_mslp_a.shape[0]
axes[2][1].set_title('f) E n='+str(ndays), loc='left')

data = neaus_mslp_a.mean(axis=0)
t, sig = stats.ttest_ind(neaus_mslp_a+lamslp_clim, laclim, axis=0, equal_var=False)
_, m = plot_mslp(data,axes[1][1])
ndays = neaus_mslp_a.shape[0]
axes[1][1].set_title('d) NE n='+str(ndays), loc='left')

data = naus_mslp_a.mean(axis=0)
t, sig = stats.ttest_ind(naus_mslp_a+lamslp_clim, laclim, axis=0, equal_var=False)
contours, m = plot_mslp(data,axes[0][1])
ndays = naus_mslp_a.shape[0]
axes[0][1].set_title('b) N n='+str(ndays), loc='left')


cax = f.add_axes([0.1,0.07,0.8,0.02])
cbar = plt.colorbar(contours,cax=cax,orientation='horizontal',ticks=np.arange(-600,601,200))
cbar.set_label('$Pa$')
f.suptitle('El Nino            La Nina', fontsize=20)
plt.savefig('mslp_hwdays_composites.eps',format='eps')
#plt.show()
print 'done'

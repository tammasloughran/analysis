import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import scipy.stats as stats
import sys

# Directory of data
directory = '/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/'
filename = directory+'EHF_heatwaves_AWAP_bp1979-2008_yearly_summer.nc'

# Load data
awapnc = nc.Dataset(filename,'r')
times = np.array(awapnc.variables['time'][:])
lats = np.array(awapnc.variables['lat'][:])
lons = np.array(awapnc.variables['lon'][:])
vrbls = ['HWA_EHF', 'HWD_EHF', 'HWF_EHF']
hwdata = np.ones((3,len(times),len(lats),len(lons)))*np.nan
for i,v in enumerate(vrbls):
    hwdata[i,...] = awapnc.variables[v][:]
hwdata[hwdata==-999.99] = np.nan

# Climatology 1979-2008
iclim = [year in range(1978,2009) for year in times]
climdata = hwdata[:,iclim,...]
clim = np.nanmean(hwdata[:,iclim,...], axis=1)

# Define El Nino and La Nina years. year is year containing December
elninoyears = [1911,1913,1914,1918,1925,1930,1941,1951,1957,1965,1969,1972,1976,1982,1987,1997,2006]
laninayears = [1909,1915,1920,1933,1942,1949,1950,1955,1970,1973,1975,1984,1988,1998,1999,2007,2010]
#elninoyears = [1972,1976,1982,1987,1997,2006]
#laninayears = [1970,1973,1975,1984,1988,1998,1999,2007,2010]
inino = [year in elninoyears for year in times]
inina = [year in laninayears for year in times]

elnino_hw = hwdata[:,inino,...]
lanina_hw = hwdata[:,inina,...]
nino_mean = np.nanmean(elnino_hw, axis=1) - clim
nina_mean = np.nanmean(lanina_hw, axis=1) - clim

masknc = nc.Dataset('/srv/ccrc/data35/z5032520/AWAP/mask/varmask.nc','r')
mask = masknc.variables['mask'][:]

levels = np.arange(-13,14,1)

fig, ax = plt.subplots(nrows=2,ncols=3,figsize=(6,4))
plt.subplots_adjust(top=0.95,bottom=0.19,left=0.075,right=0.94,hspace=0.19,wspace=0.2)

m = Basemap(ax=ax[0][0],projection='mill',
    llcrnrlon=105.,llcrnrlat=-45.,
    urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nino_mean[2], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(elnino_hw[2,...], hwdata[2,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ax[0][0].set_title('a)  Frequency',loc='left')

m = Basemap(ax=ax[1][0],projection='mill',
    llcrnrlon=105.,llcrnrlat=-45.,
    urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nina_mean[2], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(lanina_hw[2,...], hwdata[2,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[1,0,0,1],dashes=[5,700],fontsize=8)
ax[1][0].set_title('d)',loc='left')

cax = fig.add_axes([0.07,0.11,0.26,0.04])
cbar = plt.colorbar(colors,cax=cax,orientation='horizontal',ticks=[-11,-7,-3,0,3,7,11])
cbar.set_label('days')

levels = np.arange(-7,8,1)

m = Basemap(ax=ax[0][1],projection='mill',
    llcrnrlon=105.,llcrnrlat=-45.,
    urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nino_mean[1], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(elnino_hw[1,...], hwdata[1,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ax[0][1].set_title('b)  Duration',loc='left')

m = Basemap(ax=ax[1][1],projection='mill',
    llcrnrlon=105.,llcrnrlat=-45.,
    urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nina_mean[1], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(lanina_hw[1,...], hwdata[1,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,1],dashes=[5,700],fontsize=8)
ax[1][1].set_title('e)',loc='left')

cax = fig.add_axes([0.38,0.11,0.26,0.04])
cbar = plt.colorbar(colors,cax=cax,orientation='horizontal',ticks=[-6,-3,0,3,6])
cbar.set_label('days')


levels = np.arange(-28,29,4)

m = Basemap(ax=ax[0][2],projection='mill',
     llcrnrlon=105.,llcrnrlat=-45.,
     urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nino_mean[0], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(elnino_hw[0,...], hwdata[0,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
ax[0][2].set_title('c)  Amplitude',loc='left')

m = Basemap(ax=ax[1][2],projection='mill',
     llcrnrlon=105.,llcrnrlat=-45.,
     urcrnrlon=160.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = m(lns,lts)
colors = m.contourf(x,y,nina_mean[0], cmap='bwr', levels=levels)
t, sig = stats.ttest_ind(lanina_hw[0,...], hwdata[0,...], axis=0, equal_var=False, nan_policy='omit')
m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'],linewidth=0.5)
m.drawcoastlines()
m.drawparallels([-10,-20,-30,-40,-50],labels=[0,0,0,0],dashes=[5,700],fontsize=8)
m.drawmeridians([90,110,130,150,170],labels=[0,0,0,1],dashes=[5,700],fontsize=8)
ax[1][2].set_title('f)',loc='left')

cax = fig.add_axes([0.685,0.11,0.26,0.04])
cbar = plt.colorbar(colors,cax=cax,orientation='horizontal',ticks=levels[::2])
cbar.set_label('$C^{\circ2}$')


plt.savefig('awap_ENSO_HW.eps',format='eps')

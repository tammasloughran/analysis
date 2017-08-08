# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:45:42 2017

@author: Tammas Loughran
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats

# Define the ensemble members
control = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh',
           'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp',
           'vaowq','vaowr','vaows','vaowt','vaowu','vaowv','vaoww','vaowx',
           'vaowy','vaowz','vaqgi','vaqgj','vaqgk','vaowf']
elnino = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg',
          'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq',
          'vaoqr','vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy',
          'vaoqz','vaqgl','vaqgm','vaqgn','vaoqh','vaoqn']
lanina = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl',
          'vamrm','vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt',
          'vamru','vamrv','vamrw','vamrx','vamry','vamrz','vaqga','vaqgb',
          'vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']
hwdir = '/srv/ccrc/data46/z5032520/ehfheatwaves/'

modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
lsmnc = nc.Dataset(modeldir+'sftlf_fx_ACCESS1-0_historical_r0i0p0.nc', 'r')
lsm = lsmnc.variables['sftlf'][:]
lsmnc.close()
group = 'control'
print group
#Load the EF data from the control
cEF = np.empty((0,145,192))
for ens in control:
    print ens
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        qnc = nc.Dataset(ncfile,'r')
        Qe = qnc.variables['lh'][:]
        #Qh = qnc.variables['sh'][:]
        #cEF = np.append(cEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
        cEF = np.append(cEF, np.squeeze(Qe), axis=0)
        qnc.close()
cEF_clim = cEF.mean(axis=0)

group = 'elnino'
print group
#Load the EF data from the control
oEF = np.empty((0,145,192))
for ens in elnino:
    print ens
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        qnc = nc.Dataset(ncfile,'r')
        Qe = qnc.variables['lh'][:]
        #Qh = qnc.variables['sh'][:]
        oEF = np.append(oEF, np.squeeze(Qe), axis=0)
        #oEF = np.append(oEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
        qnc.close()
oEF_clim = oEF.mean(axis=0)

group = 'lanina'
print group
#Load the EF data from the control
aEF = np.empty((0,145,192))
for ens in lanina:
    print ens
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        qnc = nc.Dataset(ncfile,'r')
        Qe = qnc.variables['lh'][:]
        #Qh = qnc.variables['sh'][:]
        aEF = np.append(aEF, np.squeeze(Qe), axis=0)
        #aEF = np.append(aEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
        qnc.close()
aEF_clim = aEF.mean(axis=0)

hwnc = nc.Dataset(hwdir+'vamrd/EHF_heatwaves_ACCESS1.3_vamrd_daily.nc')
lats = hwnc.variables['lat'][:]
lons = hwnc.variables['lon'][:]
hwnc.close()

# Plot the climatology
#mp = Basemap(projection='mill',
#             llcrnrlon=110.,llcrnrlat=-48.,
#             urcrnrlon=157.,urcrnrlat=-5.)
#lns,lts = np.meshgrid(lons,lats)
#x,y = mp(lns-1,lts-1)
data = oEF_clim-aEF_clim
data = np.ma.array(data,mask=lsm<50)
#_, p = stats.ttest_ind(oEF, aEF, axis=0,equal_var=False)
#sig = p<0.05
#sig[lsm<50] = 1
#mp.contour(x,y,sig,levels=[1],colors='k')
#shade = mp.pcolormesh(x,y,data,vmin=-.3,vmax=0.3,cmap='BrBG')
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],linewidth=0)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],linewidth=0)
#mp.colorbar(shade)
#plt.show()

#f, axes = plt.subplots(nrows=2, ncols=1,figsize=(6,8))
#mp = Basemap(ax=axes[0],
#             projection='mill',
#             llcrnrlon=110.,llcrnrlat=-48.,
#             urcrnrlon=157.,urcrnrlat=-5.)
#lns,lts = np.meshgrid(lons,lats)
#x,y = mp(lns,lts)
#z = np.ma.array(cEF_clim, mask=lsm<50)
#shade = mp.pcolormesh(x,y,z,vmin=0,vmax=140,cmap='YlGnBu')
#levels = np.arange(0,141,15)
#mp.contour(x,y,z,colors='k',levels=levels,linewidths=0.3)
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],linewidth=0)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],linewidth=0)
#cbar = mp.colorbar(shade)
#cbar.set_label('$Wm^{-2}$')
#axes[0].set_title('a) Control Climatology',loc='left')
#mp = Basemap(ax=axes[1],projection='mill',
#             llcrnrlon=110.,llcrnrlat=-48.,
#             urcrnrlon=157.,urcrnrlat=-5.)
#lns,lts = np.meshgrid(lons,lats)
#x,y = mp(lns,lts)
##mp.contour(x,y,sig,levels=[1],colors='k')
#shade = mp.pcolormesh(x,y,data,vmin=-20,vmax=20,cmap='BrBG')
#levels = np.arange(-20,21,4)
#mp.contour(x,y,data,colors='k',levels=levels,linewidths=0.3)
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],linewidth=0)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],linewidth=0)
#cbar = mp.colorbar(shade)
#cbar.set_label('$Wm^{-2}$')
#axes[1].set_title('b) El Nino - La Nina',loc='left')
#plt.savefig('LHF_alldays.eps',format='eps')


#Load the metadata from the first file
hwnc = nc.Dataset(hwdir+'vamrd/EHF_heatwaves_ACCESS1.3_vamrd_daily.nc')
lats = hwnc.variables['lat'][:]
lons = hwnc.variables['lon'][:]
times = hwnc.variables['time'][:]
dates = nc.num2date(times, hwnc.variables['time'].units)
years = np.array([i.year for i in dates])
months = np.array([i.month for i in dates])
days = np.array([i.day for i in dates])
nov1 = np.where((years==2000)&(months==11)&(days==1))[0][0]
mar31 = np.where((years==2001)&(months==3)&(days==31))[0][0]
delta = int(mar31+1-nov1)

# Load the event data
event = hwnc.variables['ends'][nov1:mar31+1,...]
for ens in elnino[1:]:
    hwnc = nc.Dataset(hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_daily.nc')
    event = np.append(event, hwnc.variables['ends'][nov1:mar31+1,...],axis=0)
    hwnc.close()
event = event[:,:,(lons>=112.5)&(lons<=155.625)]
event = event[:,(lats>=-43.75)&(lats<=-10),:]

# trim lons and lats for australia
lats2 = lats[(lats>=-43.75)&(lats<=-10)]
lons2 = lons[(lons>=112.5)&(lons<=155.625)]

# define region of interest
def region_mask(lon,lat,size=7.5):
    region = np.zeros(event.shape[1:])
    for i,x in enumerate(lons2):
        for j,y in enumerate(lats2):
            if (x>=lon)&(x<=lon+size)&(y<=lat)&(y>=lat-size):
                region[j,i] = 1
    return np.repeat(region[None,...], event.shape[0], axis=0)


#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)
seaus = region_mask(141.,-31.)
neaus = region_mask(139.,-18.)
naus = region_mask(129.,-12)
eaus = region_mask(145.5,-24)


def plot_ef(data,ax,ndays=0):
    data = np.ma.array(data, mask=lsm<50)
    m = Basemap(ax=ax,projection='mill',
                llcrnrlon=110.,llcrnrlat=-48.,
                urcrnrlon=157.,urcrnrlat=-5.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    cont = m.pcolormesh(x,y,data,cmap='PuOr',vmin=-50,vmax=50)
    m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'])
    levels = np.arange(-50,51,10)
    cnt = m.contour(x,y,data,colors='k',linewidths=0.3,levels=levels)
    for c in cnt.collections:
        if c.get_linestyle() == [(None, None)]:
            continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    m.drawcoastlines()
    return cont

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(5.5,7.75),sharex=True,sharey=True)

# create index using the region of interest
event.mask = np.logical_not(seaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[3][0])
ndays = index.sum()
axes[3][0].set_title('g) n='+str(ndays), loc='left')

event.mask = np.logical_not(eaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[2][0])
ndays = index.sum()
axes[2][0].set_title('e) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[1][0])
ndays = index.sum()
axes[1][0].set_title('c) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[0][0])
ndays = index.sum()
axes[0][0].set_title('a) n='+str(ndays), loc='left')

#Load the metadata from the first file
hwnc = nc.Dataset(hwdir+'vamre/EHF_heatwaves_ACCESS1.3_vamre_daily.nc')
lats = hwnc.variables['lat'][:]
lons = hwnc.variables['lon'][:]
times = hwnc.variables['time'][:]
dates = nc.num2date(times, hwnc.variables['time'].units)
years = np.array([i.year for i in dates])
months = np.array([i.month for i in dates])
days = np.array([i.day for i in dates])
nov1 = np.where((years==2000)&(months==11)&(days==1))[0][0]
mar31 = np.where((years==2001)&(months==3)&(days==31))[0][0]
delta = int(mar31+1-nov1)

# Load the event data
event = hwnc.variables['event'][nov1:mar31+1,...]
for i,ens in enumerate(lanina[1:]):
    hwnc = nc.Dataset(hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_daily.nc')
    event = np.append(event, hwnc.variables['event'][nov1:mar31+1,...],axis=0)
    hwnc.close()
event = event[:,:,(lons>=112.5)&(lons<=155.625)]
event = event[:,(lats>=-43.75)&(lats<=-10),:]

# trim lons and lats for australia
lats2 = lats[(lats>=-43.75)&(lats<=-10)]
lons2 = lons[(lons>=112.5)&(lons<=155.625)]

event.mask = np.logical_not(seaus)
index = event.sum(axis=1).sum(axis=1)
index = index>10
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[3][1])
ndays = index.sum()
axes[3][1].set_title('h) n='+str(ndays), loc='left')
    
event.mask = np.logical_not(eaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[2][1])
ndays = index.sum()
axes[2][1].set_title('f) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
plot_ef(data,axes[1][1])
ndays = index.sum()
axes[1][1].set_title('d) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = 1
mesh = plot_ef(data,axes[0][1])
ndays = index.sum()
axes[0][1].set_title('b) n='+str(ndays), loc='left')

cax = f.add_axes([0.1,0.07,0.8,0.02])
cbar = plt.colorbar(mesh,cax=cax,orientation='horizontal',ticks=np.arange(-50,51,10))
cbar.set_label('$Wm^{-2}$')
f.suptitle('El Nino            La Nina', fontsize=20)
plt.savefig('LHF_composites.eps',format='eps')
print 'done'
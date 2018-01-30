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
import matplotlib
matplotlib.rcParams['contour.negative_linestyle']= 'dashed'
import sys

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
#cEF = np.empty((0,145,192))
#cQe = np.empty((0,145,192))
#cQh = np.empty((0,145,192))
#for ens in control:
#    print ens
#    for month in ['11','12','01','02','03']:
#        if month=='11' or month=='12': 
#            year = '2000'
#        else:
#            year = '2001'
#        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
#        qnc = nc.Dataset(ncfile,'r')
#        Qe = qnc.variables['lh'][:]
#        Qh = qnc.variables['sh'][:]
#        cEF = np.append(cEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
#        cQe = np.append(cQe, np.squeeze(Qe), axis=0)
#        cQh = np.append(cQh, np.squeeze(Qh), axis=0)
#        qnc.close()
cEF = np.load('cEF.npy')
cEF_clim = cEF.mean(axis=0)
cQe = np.load('cQe.npy')
cQe_clim = cQe.mean(axis=0)
cQh = np.load('cQh.npy')
cQh_clim = cQh.mean(axis=0)

group = 'elnino'
print group
#Load the EF data from the control
#oEF = np.empty((0,145,192))
#oQe = np.empty((0,145,192))
#oQh = np.empty((0,145,192))
#for ens in elnino:
#    print ens
#    for month in ['11','12','01','02','03']:
#        if month=='11' or month=='12': 
#            year = '2000'
#        else:
#            year = '2001'
#        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
#        qnc = nc.Dataset(ncfile,'r')
#        Qe = qnc.variables['lh'][:]
#        Qh = qnc.variables['sh'][:]
#        oEF = np.append(oEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
#        oQe = np.append(oQe, np.squeeze(Qe), axis=0)
#        oQh = np.append(oQh, np.squeeze(Qh), axis=0)
#        qnc.close()
oEF = np.load('oEF.npy')
oEF_clim = oEF.mean(axis=0)
oQe = np.load('oQe.npy')
oQe_clim = oQe.mean(axis=0)
oQh = np.load('oQh.npy')
oQh_clim = oQh.mean(axis=0)

group = 'lanina'
print group
#Load the EF data from the control
#aEF = np.empty((0,145,192))
#aQe = np.empty((0,145,192))
#aQh = np.empty((0,145,192))
#for ens in lanina:
#    print ens
#    for month in ['11','12','01','02','03']:
#        if month=='11' or month=='12': 
#            year = '2000'
#        else:
#            year = '2001'
#        ncfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
#        qnc = nc.Dataset(ncfile,'r')
#        Qe = qnc.variables['lh'][:]
#        Qh = qnc.variables['sh'][:]
#        aEF = np.append(aEF, np.squeeze(Qe/(Qe+Qh)), axis=0)
#        aQe = np.append(aQe, np.squeeze(Qe), axis=0)
#        aQh = np.append(aQh, np.squeeze(Qh), axis=0)
#        qnc.close()
aEF = np.load('aEF.npy')
aEF_clim = aEF.mean(axis=0)
aQe = np.load('aQe.npy')
aQe_clim = aQe.mean(axis=0)
aQh = np.load('aQh.npy')
aQh_clim = aQh.mean(axis=0)



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
EFdiff = oEF_clim-aEF_clim
EFdiff[lsm<50] = 0
Qediff = oQe_clim-aQe_clim
Qediff = np.ma.array(Qediff,mask=lsm<50)
Qhdiff = oQh_clim-aQh_clim
Qhdiff = np.ma.array(Qhdiff,mask=lsm<50)
#_, p = stats.ttest_ind(oEF, aEF, axis=0,equal_var=False)
#sig = p<0.05
#sig[lsm<50] = 1
#mp.contour(x,y,sig,levels=[1],colors='k')
#shade = mp.pcolormesh(x,y,data,vmin=-.3,vmax=0.3,cmap='BrBG')
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],linestyle='none')
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],linestyle='none')
#mp.colorbar(shade)
#plt.show()

f, axes = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
## The Evaporative fraction clim
#mp = Basemap(ax=axes[0,0],projection='mill',
#             llcrnrlon=110.,llcrnrlat=-48.,
#             urcrnrlon=157.,urcrnrlat=-5.)
#lns,lts = np.meshgrid(lons,lats)
#x,y = mp(lns-1,lts-1)
#xc,yc = mp(lns,lts)
#z = np.ma.array(cEF_clim, mask=lsm<50)
#shade = mp.pcolormesh(x,y,z,vmin=0,vmax=1,cmap='terrain_r')
#levels = np.arange(0,1,0.1)
#mp.contour(xc,yc,z,colors='k',levels=levels,linewidths=0.6)
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
#mp.colorbar(shade)
#axes[0,0].set_title('a) EF Control Climatology',loc='left')
# Evaporative fraction difference
#mp = Basemap(ax=axes[1,0],projection='mill',
#             llcrnrlon=110.,llcrnrlat=-48.,
#             urcrnrlon=157.,urcrnrlat=-5.)
#lns,lts = np.meshgrid(lons,lats)
#x,y = mp(lns-1,lts-1)
#xc,yc = mp(lns,lts)
#mp.contour(x,y,sig,levels=[1],colors='k')
#shade = mp.pcolormesh(x,y,EFdiff,vmin=-.3,vmax=0.3,cmap='BrBG')
#levels = np.arange(-0.3,0.35,0.06)
#mp.contour(xc,yc,EFdiff,colors='k',levels=levels,linewidths=0.6)
#mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
#cbar = mp.colorbar(shade)
#axes[1,0].set_title('d) EF El Nino - La Nina',loc='left')
# Sensible heat clim
mp = Basemap(ax=axes[0,0],projection='cyl',
             llcrnrlon=110.,llcrnrlat=-48.,
             urcrnrlon=157.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = mp(lns-1,lts-1)
xc,yc = mp(lns,lts)
zm = cQh_clim
z = np.ma.array(cQh_clim, mask=lsm<50)
shade = mp.pcolormesh(x,y,z,vmin=0,vmax=135,cmap='YlOrRd')
levels = np.arange(0,135,15)
mp.contour(xc,yc,z,colors='k',levels=levels,linewidths=0.6)
mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
cbar = mp.colorbar(shade)
cbar.set_label('$Wm^{-2}$')
axes[0,0].set_title('a) Q$_{H}$ Control Climatology',loc='left')
# Sensible heat diff
mp = Basemap(ax=axes[1,0],projection='cyl',
             llcrnrlon=110.,llcrnrlat=-48.,
             urcrnrlon=157.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = mp(lns-1,lts-1)
xc,yc = mp(lns,lts)
z = np.ma.array(Qhdiff, mask=lsm<50)
shade = mp.pcolormesh(x,y,z,vmin=-30,vmax=30,cmap='bwr')
levels = np.arange(-30,30,5)
mp.contour(xc,yc,z,colors='k',levels=levels,linewidths=0.6)
mp.drawcoastlines()
mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
cbar = mp.colorbar(shade)
cbar.set_label('$Wm^{-2}$')
axes[1,0].set_title('c) Q$_{H}$ El Nino - La Nina',loc='left')
# latent heat clim
mp = Basemap(ax=axes[0,1],projection='cyl',
             llcrnrlon=110.,llcrnrlat=-48.,
             urcrnrlon=157.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = mp(lns-1,lts-1)
xc,yc = mp(lns,lts)
z = np.ma.array(cQe_clim, mask=lsm<50)
shade = mp.pcolormesh(x,y,z,vmin=0,vmax=135,cmap='YlGnBu')
levels = np.arange(0,135,15)
mp.contour(xc,yc,z,colors='k',levels=levels,linewidths=0.6)
mp.drawcoastlines()
#mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
cbar = mp.colorbar(shade)
cbar.set_label('$Wm^{-2}$')
axes[0,1].set_title('b) Q$_{E}$ Control Climatology',loc='left')
# Latent heat diff
mp = Basemap(ax=axes[1,1],projection='cyl',
             llcrnrlon=110.,llcrnrlat=-48.,
             urcrnrlon=157.,urcrnrlat=-5.)
lns,lts = np.meshgrid(lons,lats)
x,y = mp(lns-1,lts-1)
xc,yc = mp(lns,lts)
z = np.ma.array(Qediff, mask=lsm<50)
shade = mp.pcolormesh(x,y,z,vmin=-20,vmax=20,cmap='BrBG')
levels = np.arange(-20,20,4)
mp.contour(xc,yc,z,colors='k',levels=levels,linewidths=0.6)
#t, sig = stats.ttest_ind(oQe, aQe, axis=0, equal_var=False)
#mp.contourf(xc,yc,sig<0.05, 1, colors='none', hatches=[None,'xxx'])
mp.drawcoastlines()
mp.drawmeridians(np.arange(110,151,10),labels=[0,0,0,1],dashes=[5,700],fontsize=9)
#mp.drawparallels(np.arange(-40,-5,10),labels=[1,0,0,0],dashes=[5,700],fontsize=9)
cbar = mp.colorbar(shade)
cbar.set_label('$Wm^{-2}$')
axes[1,1].set_title('d) Q$_{E}$ El Nino - La Nina',loc='left')
plt.subplots_adjust(hspace=0.01,wspace=0.25)
#plt.show()
plt.savefig('EF_alldays.eps',format='eps')
print 'figure done'
sys.exit()

#Load the metadata from the first file
hwnc = nc.Dataset(hwdir+'vamrd/EHF_heatwaves_ACCESS1.3_vamrd_daily.nc')
lats = hwnc.variables['lat'][:]
lons = hwnc.variables['lon'][:]
times = hwnc.variables['time'][:]
dates = nc.num2date(times, hwnc.variables['time'].units)
years = np.array([i.year for i in dates])
months = np.array([i.month for i in dates])
days = np.array([i.day for i in dates])
nov1 = np.where((years==2000)&(months==11)&(days==1))[0]
mar31 = np.where((years==2001)&(months==3)&(days==31))[0]
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
    #cont = m.pcolormesh(x,y,data,cmap='terrain_r',vmin=0,vmax=1)
    cont = m.pcolormesh(x,y,data,cmap='BrBG',vmin=-.4,vmax=.4)
    #levels = np.arange(0,1,.1)
    m.contourf(x,y,sig<0.05, 1, colors='none', hatches=[None,'xxx'])
    levels = np.arange(-.4,.48,0.08)
    cnt = m.contour(x,y,data,colors='k',linewidths=0.3,levels=levels)
    for c in cnt.collections:
        if c.get_linestyle() == [(None, None)]:
            continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    m.drawcoastlines()
    #m.drawparallels(lats,linestyle='none',labels=[1,0,0,1])
    #m.drawmeridians(lons,linestyle='none',labels=[1,0,0,1])
    return cont

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(5.5,7.75),sharex=True,sharey=True)

# create index using the region of interest
event.mask = np.logical_not(seaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#only = oEF[index,...][:176]
abs_data = oEF[index,...] 
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[3][0])
ndays=index.sum()
#ndays = 176
axes[3][0].set_title('g) n='+str(ndays), loc='left')

event.mask = np.logical_not(eaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#only = oEF[index,...][:95]
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[2][0])
ndays=index.sum()
#ndays = 95
axes[2][0].set_title('e) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#only = oEF[index,...][:87]
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[1][0])
ndays=index.sum()
#ndays = 87
axes[1][0].set_title('c) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#only = oEF[index,...][:50]
abs_data = oEF[index,...]
data = oEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[0][0])
ndays=index.sum()
#ndays = 50
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
nov1 = np.where((years==2000)&(months==11)&(days==1))[0]
mar31 = np.where((years==2001)&(months==3)&(days==31))[0]
delta = int(mar31+1-nov1)

abs_data = 0

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
sig[lsm<50] = np.nan
plot_ef(data,axes[3][1])
ndays=index.sum()
axes[3][1].set_title('h) n='+str(ndays), loc='left')
    
event.mask = np.logical_not(eaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[2][1])
ndays=index.sum()
axes[2][1].set_title('f) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
plot_ef(data,axes[1][1])
ndays=index.sum()
axes[1][1].set_title('d) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>9
abs_data = aEF[index,...]
data = aEF[index,...].mean(axis=0) - cEF_clim
data[lsm<50] = 0
t, sig = stats.ttest_ind(abs_data, cEF, axis=0, equal_var=False)
sig[lsm<50] = np.nan
mesh = plot_ef(data,axes[0][1])
ndays=index.sum()
axes[0][1].set_title('b) n='+str(ndays), loc='left')

cax = f.add_axes([0.1,0.07,0.8,0.02])
plt.colorbar(mesh,cax=cax,orientation='horizontal')
f.suptitle('El Nino            La Nina', fontsize=20)
plt.savefig('EF_composites.eps',format='eps')

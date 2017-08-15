# -*- coding: utf-8 -*-
"""
Created on Wed May  3 14:42:49 2017

@author: Tammas Loughran
"""


import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats

# define the ensembles
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

# Load the pressure climatology
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'control'
pr = np.empty((0,145,192))
for ens in control:
    print ens
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        prnc = nc.Dataset(prfile,'r')
        pr = np.append(pr, np.squeeze(prnc.variables['p_1'][:]).mean(axis=0)[None,...], axis=0)
        prnc.close()
prc = pr
pclim = np.ones((5,)+pr.shape[1:])*np.nan
for i in xrange(5):
    pclim[i,...] = pr[i::5].mean(axis=0)
    
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

# Load the pressure data
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'elnino'
pr = np.empty((0,145,192))
for ens in elnino:
    print ens
    for i,month in enumerate(['11','12','01','02','03']):
        prclim = pclim[i,...]
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        prnc = nc.Dataset(prfile,'r')
        pr = np.append(pr, np.squeeze(prnc.variables['p_1'][:]-pclim[i]), axis=0)
        prnc.close()


def plot_pr(data,p,ax,ndays=0):
    m = Basemap(ax=ax,projection='mill',
            llcrnrlon=70.,llcrnrlat=-60.,
            urcrnrlon=220.,urcrnrlat=20.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    cont = m.contourf(x,y,data,cmap='bwr',levels=range(-700,800,100),extend='both')
    m.contourf(x,y,p<0.03,1,colors='none',hatches=[None,'xx'])
    m.drawcoastlines()
    #m.drawparallels(lats,labels=[1,0,0,1])
    #m.drawmeridians(lons,labels=[1,0,0,1])
    return cont
    
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

f, axes = plt.subplots(nrows=4, ncols=2,figsize=(5.5,7.75),sharex=True,sharey=True)

# create index using the region of interest
event.mask = np.logical_not(seaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[3][0])
ndays=index.sum()
axes[3][0].set_title('g) n='+str(ndays), loc='left')

event.mask = np.logical_not(eaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[2][0])
ndays=index.sum()
axes[2][0].set_title('e) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[1][0])
ndays=index.sum()
axes[1][0].set_title('c) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[0][0])
ndays=index.sum()
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

# Load the pressure data
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'lanina'
pr = np.empty((0,145,192))
for ens in lanina:
    print ens
    for i,month in enumerate(['11','12','01','02','03']):
        prclim = pclim[i,...]
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        prnc = nc.Dataset(prfile,'r')
        pr = np.append(pr, np.squeeze(prnc.variables['p_1'][:]-pclim[i]), axis=0)
        prnc.close()

# define region of interest
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)
seaus = region_mask(141.,-31.)
neaus = region_mask(139.,-18.)
naus = region_mask(129.,-12)
eaus = region_mask(145.5,-24)

# create index using the region of interest (& tasman pr>1020hpa)
event.mask = np.logical_not(seaus)
index = event.sum(axis=1).sum(axis=1)
index = index>10
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[3][1])
ndays=index.sum()
axes[3][1].set_title('h) n='+str(ndays), loc='left')
    
event.mask = np.logical_not(eaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[2][1])
ndays=index.sum()
axes[2][1].set_title('f) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
plot_pr(pr[index,...].mean(axis=0),p,axes[1][1])
ndays=index.sum()
axes[1][1].set_title('d) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>9
_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
cont = plot_pr(pr[index,...].mean(axis=0),p,axes[0][1])
ndays=index.sum()
axes[0][1].set_title('b) n='+str(ndays), loc='left')

cax = f.add_axes([0.1,0.07,0.8,0.02])
plt.colorbar(cont,cax=cax,orientation='horizontal')
cax.set_xlabel('Pa')
f.suptitle('El Nino            La Nina', fontsize=20)
plt.show()
#plt.savefig('mslp_composites.eps',format='eps')

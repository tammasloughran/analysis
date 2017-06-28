# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:51:54 2017

@author: Tammas Loughran
"""

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
import pdb

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
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        prnc = nc.Dataset(prfile,'r')
        pr = np.append(pr, np.squeeze(prnc.variables['p_1'][:]), axis=0)
        prnc.close()
pclim = pr.mean(axis=0)

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
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        prnc = nc.Dataset(prfile,'r')
        pr = np.append(pr, np.squeeze(prnc.variables['p_1'][:]), axis=0)
        prnc.close()
pclim = pr.mean(axis=0)

def plot_pr(data,ndays=0):
    m = Basemap(projection='mill',
            llcrnrlon=70.,llcrnrlat=-60.,
            urcrnrlon=220.,urcrnrlat=20.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    m.contourf(x,y,data,cmap='bwr',levels=range(-700,800,100),extend='both')
    m.drawcoastlines()
    plt.title('n='+str(ndays))
    cb = m.colorbar()
    m.drawparallels(lats,linewidth=0)
    m.drawmeridians(lons,linewidth=0)
    plt.show()
    
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

# create index using the region of interest (& tasman pr>1020hpa)
event.mask = np.logical_not(seaus)
index = event.sum(axis=1).sum(axis=1)
index = index>10
plot_pr(pr[index,...].mean(axis=0)-pclim,ndays=index.sum())
    
event.mask = np.logical_not(eaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
plot_pr(pr[index,...].mean(axis=0)-pclim,ndays=index.sum())

event.mask = np.logical_not(neaus)
index = event.sum(axis=1).sum(axis=1)
index = index>9
plot_pr(pr[index,...].mean(axis=0)-pclim,ndays=index.sum())

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>9
plot_pr(pr[index,...].mean(axis=0)-pclim,ndays=index.sum())
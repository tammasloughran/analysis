# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:54:06 2017

@author: tammas
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
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

# Load the wind climatology
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'control'
u = np.empty((0,144,192))
v = np.empty((0,144,192))
for ens in control:
    print ens
    for month in ['11','12','01','02','03']:
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        uvfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        uvnc = nc.Dataset(uvfile,'r')
        u = np.append(u, np.squeeze(uvnc.variables['u'][:]).mean(axis=0)[None,...], axis=0)
        v = np.append(u, np.squeeze(uvnc.variables['v'][:]).mean(axis=0)[None,...], axis=0)
uc = u
vc = v
ucclim = np.ones((5,)+u.shape[1:])*np.nan
vcclim = np.ones((5,)+v.shape[1:])*np.nan
for i in xrange(5):
    ucclim[i,...] = uc[i::5].mean(axis=0)
    vcclim[i,...] = vc[i::5].mean(axis=0)



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
#lats1 = lats1[(lats>=-43.75)&(lats<=-10)]
#lons1 = lons1[(lons>=112.5)&(lons<=155.625)]
lats = uvnc.variables['latitude_1'][:]
lons = uvnc.variables['longitude_1'][:]

# Load the wind data
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'elnino'
u = np.empty((0,144,192))
v = np.empty((0,144,192))
for ens in elnino:
    print ens
    for i,month in enumerate(['11','12','01','02','03']):
        uclim = ucclim[i,...]
        vclim = vcclim[i,...]
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        uvfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        uvnc = nc.Dataset(uvfile,'r')
        u = np.append(u, np.squeeze(uvnc.variables['u'][:]-uclim[i]), axis=0)
        v = np.append(v, np.squeeze(uvnc.variables['v'][:]-vclim[i]), axis=0)
        uvnc.close()


def plot_pr(u,v,ax,ndays=0):
    m = Basemap(ax=ax,projection='merc',
            llcrnrlon=100.,llcrnrlat=-60.,
            urcrnrlon=200.,urcrnrlat=20.)
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    ugrid,newlons = shiftgrid(180.,ucomp,lons,start=False)
    vgrid,newlons = shiftgrid(180.,vcomp,lons,start=False)
    uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,lats,20,20,returnxy=True)
    Q = m.quiver(xx,yy,uproj,vproj,scale=100)
    #m.contourf(x,y,p<0.03,1,colors='none',hatches=[None,'xx'])
    m.drawcoastlines()
    #m.drawparallels(lats,labels=[1,0,0,1])
    #m.drawmeridians(lons,labels=[1,0,0,1])
    return Q
    
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
f.suptitle('El Nino            La Nina', fontsize=20)
# create index using the region of interest
event.mask = np.logical_not(seaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
quivers = plot_pr(ucomp,vcomp,axes[3][0])
ndays=index.sum()
axes[3][0].set_title('g) n='+str(ndays), loc='left')
plt.quiverkey(quivers,0.1,-0.1,10,'10 m/s',labelpos='E')

event.mask = np.logical_not(eaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[2][0])
ndays=index.sum()
axes[2][0].set_title('e) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[1][0])
ndays=index.sum()
axes[1][0].set_title('c) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[0][0])
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

# Load the wind data
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
group = 'lanina'
u = np.empty((0,144,192))
v = np.empty((0,144,192))
for ens in lanina:
    print ens
    for i,month in enumerate(['11','12','01','02','03']):
        uclim = ucclim[i,...]
        vclim = vcclim[i,...]
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        uvfile = modeldir+group+'/'+ens+'/'+ens+'a.pe'+year+'-'+month+'.nc'
        uvnc = nc.Dataset(uvfile,'r')
        u = np.append(u, np.squeeze(uvnc.variables['u'][:]-uclim[i]), axis=0)
        v = np.append(v, np.squeeze(uvnc.variables['v'][:]-vclim[i]), axis=0)
        uvnc.close()

# define region of interest
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)
seaus = region_mask(141.,-31.)
neaus = region_mask(139.,-18.)
naus = region_mask(129.,-12)
eaus = region_mask(145.5,-24)

# create index using the region of interest (& tasman pr>1020hpa)
event.mask = np.logical_not(seaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
quivers = plot_pr(ucomp,vcomp,axes[3][1])
ndays=index.sum()
axes[3][1].set_title('g) n='+str(ndays), loc='left')
plt.quiverkey(quivers,0.1,-0.1,10,'10 m/s',labelpos='E')

event.mask = np.logical_not(eaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[2][1])
ndays=index.sum()
axes[2][1].set_title('e) n='+str(ndays), loc='left')

event.mask = np.logical_not(neaus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[1][1])
ndays=index.sum()
axes[1][1].set_title('c) n='+str(ndays), loc='left')

event.mask = np.logical_not(naus)
event.mask[event<0] = 1
index = event.sum(axis=1).sum(axis=1)
index = index>10
#_, p = stats.ttest_ind(pr[index,...], pr, axis=0, equal_var=False)
ucomp = u[index,...].mean(axis=0)
vcomp = v[index,...].mean(axis=0)
plot_pr(ucomp,vcomp,axes[0][1])
ndays=index.sum()
axes[0][1].set_title('a) n='+str(ndays), loc='left')

plt.savefig('winds_composites.eps',format='eps')
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:47:36 2017

@author: Tammas Loughran

Investigation of El Nino ensemble members that have very few heatwaves in 
the northeast of asutralia.
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys
from pyclimate import svdeofs
import scipy.stats as stats
nino_few = ['vaoqa','vaoqd','vaoqe','vaoqf','vaoqg','vaoql','vaoqm','vaoqn',
            'vaoqt','vaoqu','vaoqv']
nino_many = ['vamrd','vaoqb','vaoqc','vaoqh','vaoqi','vaoqj','vaoqk','vaoqo',
             'vaoqp','vaoqq','vaoqr','vaoqs','vaoqw','vaoqx','vaoqy','vaoqz',
             'vaqgl','vaqgm','vaqgn']
control = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh',
           'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp',
           'vaowq','vaowr','vaows','vaowt','vaowu','vaowv','vaoww','vaowx',
           'vaowy','vaowz','vaqgi','vaqgj','vaqgk','vaowf']
lanina = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl',
          'vamrm','vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt',
          'vamru','vamrv','vamrw','vamrx','vamry','vamrz','vaqga','vaqgb',
          'vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']

# Load the MSLP data
#spindir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/vamrb/'
#monthlync = nc.MFDataset(spindir+'vamrba.pa*')
#lons = monthlync.variables['longitude'][:]
#lats = monthlync.variables['latitude'][:]
#tunits = monthlync.variables['t'].units
#times = monthlync.variables['t'][:]
#dates = nc.num2date(times,units=tunits)
#mslp_spin = monthlync.variables['p'][-30*12:,0,...]
#monthlync.close()
#
## Calculate and remove the seasonal cycle
#seasc = np.ones((12,)+mslp_spin.shape[1:])*np.nan
#for i in xrange(12):
#    seasc[i,...] = mslp_spin[i::12].mean(axis=0)
#for i in xrange(30*12):
#    mslp_spin[i,...] = mslp_spin[i,...] - seasc[np.mod(i,12),...]
#
## Weight by the square root of the cosine of the latitude
#x,y = np.meshgrid(lons,lats)
#coslat = np.sqrt(np.cos(y*np.pi/180.))
#for i in xrange(30*12):
#    mslp_spin[i,...] = mslp_spin[i,...]*coslat
#
## Cut out the SAM region
#mslp_spin = mslp_spin[:,(lats<=-25)&(lats>-75),:]
#
## Zonal Average
#mslp_spin = mslp_spin.mean(axis=2)
#
## Eofs of mslp
#svdc = svdeofs.SVDEOFs(mslp_spin)
#eofs = svdc.eofs(pcscaling=1)
#eofs_unsc = svdc.eofs(pcscaling=0)
#pcs = svdc.pcs(pcscaling=1)
#eigen = svdc.eigenvalues()
#sam = pcs[:,0]*-1
#pattern = eofs[:,0]*-1
#
## SAM pattern
#plt.plot(lats[(lats<=-25)&(lats>-75)],pattern)
#plt.xlabel('latitude')
#plt.ylabel('hPa')
#plt.show()
#
## SAM index
#plt.plot(sam)
#plt.xlabel('time')
#plt.ylabel('SAM')
#plt.show()
#
## Repeat with ensemble members projected on to the unscaled eofs.
#sam_nino = np.ones((30,24))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/'
#for j,ens in enumerate(nino_few+nino_many):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    mslp_ens = np.squeeze(ncfile.variables['p'][:])
#    for i in xrange(24):
#        mslp_ens[i,...] = mslp_ens[i,...] - seasc[np.mod(i,12),...]
#    mslp_ens = mslp_ens*coslat
#    mslp_ens = mslp_ens[:,(lats<=-25)&(lats>-75),:]
#    mslp_ens = mslp_ens.mean(axis=2)
#    pcs_ens = np.dot(mslp_ens,eofs) #[T,S]x[S,N]
#    pcs_ens = pcs_ens/pcs_ens.std(axis=0)
#    sam_nino[j,:] = pcs_ens[:,0]
#    print ens, ' DJF SAM: ', sam_nino[j,11:14].mean()
#print 'Few heatwaves group DJF SAM: ', sam_nino[:11,11:14].mean()
#print 'Many heatwaves group DJF SAM: ', sam_nino[11:,11:14].mean()
#print 'Nino DJF SAM mean: ', sam_nino[:,11:14].mean()
#djfi = np.array([1,1,0,0,0,0,0,0,0,0,0,1])
#for i in xrange(29): djfi = np.concatenate((djfi,np.array([1,1,0,0,0,0,0,0,0,0,0,1])))
#print 'Control DJF SAM mean:', sam[djfi.astype('bool')].mean()
#
#sam_nina = np.ones((30,24))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/'
#for j,ens in enumerate(lanina):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    mslp_ens = np.squeeze(ncfile.variables['p'][:])
#    for i in xrange(24):
#        mslp_ens[i,...] = mslp_ens[i,...] - seasc[np.mod(i,12),...]
#    mslp_ens = mslp_ens*coslat
#    mslp_ens = mslp_ens[:,(lats<=-25)&(lats>-75),:]
#    mslp_ens = mslp_ens.mean(axis=2)
#    pcs_ens = np.dot(mslp_ens,eofs)
#    pcs_ens = pcs_ens/pcs_ens.std(axis=0)
#    sam_nina[j,:] = pcs_ens[:,0]
#    print ens, ' DJF SAM: ', sam_nina[j,11:14].mean()
#print 'Nina DJF SAM mean: ', sam_nina[:,11:14].mean()
#
## Validate using the pressure difference index AOI.
#print "Repeat using the difference based SAM index of Gong & Wang (1999)"
#sam_nino = np.ones((30,24))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/'
#for j,ens in enumerate(nino_few+nino_many):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    mslp_ens = np.squeeze(ncfile.variables['p'][:])
#    mslp_ens = mslp_ens.mean(axis=2)
#    mslp_ens = (mslp_ens-mslp_ens.mean(axis=0))/mslp_ens.std(axis=0)
#    p40 = mslp_ens[:,np.where(lats==-40)[0][0]]
#    p65 = mslp_ens[:,np.where(lats==-65)[0][0]]
#    sam_nino[j,:] = p40-p65
#    print ens, ' DJF SAM: ', sam_nino[j,11:14].mean()
#print 'Few heatwaves group DJF SAM: ', sam_nino[:11,11:14].mean()
#print 'Many heatwaves group DJF SAM: ', sam_nino[11:,11:14].mean()
#print 'Nino DJF SAM mean: ', sam_nino[:,11:14].mean()
#
#sam_nina = np.ones((30,24))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/'
#for j,ens in enumerate(lanina):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    mslp_ens = np.squeeze(ncfile.variables['p'][:])
#    mslp_ens = mslp_ens.mean(axis=2)
#    mslp_ens = (mslp_ens-mslp_ens.mean(axis=0))/mslp_ens.std(axis=0)
#    p40 = mslp_ens[:,np.where(lats==-40)[0][0]]
#    p65 = mslp_ens[:,np.where(lats==-65)[0][0]]
#    sam_nina[j,:] = p40-p65
#    print ens, ' DJF SAM: ', sam_nina[j,11:14].mean()
#print 'Nina DJF SAM mean: ', sam_nina[:,11:14].mean()

################################
# Soil mousture & surface fluxes
################################

spindir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/vamrb/'
monthlync = nc.MFDataset(spindir+'vamrba.pa*')
lons = monthlync.variables['longitude'][:]
lats = monthlync.variables['latitude'][:]
tunits = monthlync.variables['t'].units
times = monthlync.variables['t'][:]
dates = nc.num2date(times,units=tunits)
monthlync.close()
djfi = np.array([1,1,0,0,0,0,0,0,0,0,0,1])
for i in xrange(29): djfi = np.concatenate((djfi,np.array([1,1,0,0,0,0,0,0,0,0,0,1])))
mlons = lons[(lons<155)&(lons>110)] - 1
mlats = lats[(lats<-7)&(lats>-47)] - 1
x,y = np.meshgrid(mlons,mlats)

msknc = nc.Dataset('/srv/ccrc/data46/z5032520/modelout/ACCESS/sftlf_fx_ACCESS1-0_historical_r0i0p0.nc','r')
mask = msknc.variables['sftlf'][...]
mask = mask[:,(lons<155)&(lons>110)]
mask = mask[(lats<-7)&(lats>-47),:]

# Figure for sensible, latent heat flux and surface soil moisture.
fig, axes = plt.subplots(nrows=1,ncols=3, figsize=(11,5))

sh_nino = np.ones((30,len(mlats),len(mlons)))*np.nan
ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/'
for j,ens in enumerate(nino_few+nino_many):
    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
    sh_ens = np.squeeze(ncfile.variables['sh'][:,0,...])
    sh_ens = sh_ens[11:13,...].mean(axis=0)
    sh_ens = sh_ens[(lats<-7)&(lats>-47),:]
    sh_ens = sh_ens[:,(lons<155)&(lons>110)]
    sh_nino[j,...] = sh_ens

#sh_nina = np.ones((30,len(mlats),len(mlons)))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/'
#for j,ens in enumerate(lanina):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    sh_ens = np.squeeze(ncfile.variables['sh'][:,0,...])
#    sh_ens = sh_ens[11:13,...].mean(axis=0)
#    sh_ens = sh_ens[(lats<-7)&(lats>-47),:]
#    sh_ens = sh_ens[:,(lons<155)&(lons>110)]
#    sh_nina[j,...] = sh_ens
mp = Basemap(ax=axes[0],projection='cyl',llcrnrlon=110, llcrnrlat=-47, urcrnrlon=155, urcrnrlat=-7)
mp.drawmeridians([120,130,140,150],labels=[1,0,0,1],dashes=[5,700],fontsize=9)
mp.drawparallels([-10,-20,-30,-40],lables=[1,0,0,1],dashes=[5,700],fontsize=9)
xx,yy = mp(x,y)
dat = sh_nino[:11].mean(axis=0)-sh_nino[11:].mean(axis=0)
dat[mask==0] = 0
dat = np.ma.array(dat, mask=mask<50)
mesh = mp.pcolormesh(xx,yy, dat, cmap='bwr',vmin=-15,vmax=15)
levels = np.arange(-15,16,5)
cont = mp.contour(xx,yy,dat,levels=levels,colors='k',linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
#_, p = stats.ttest_ind(sh_nino[:11], sh_nino[11:], axis=0, equal_var=False)
#mp.contourf(xx,yy,p<0.05,1,colors='none',hatches=[None,'xx'])
cbar = mp.colorbar(mesh,location='bottom',ticks=levels,pad=0.2)
cbar.set_label('Wm$^{-2}$')
mp.drawcoastlines()
plt.sca(axes[0])
plt.title('a) Q$_{H}$',loc='left')

#-------------------------------------
lh_nino = np.ones((30,len(mlats),len(mlons)))*np.nan
ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/'
for j,ens in enumerate(nino_few+nino_many):
    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
    lh_ens = np.squeeze(ncfile.variables['lh'][:,0,...])
    lh_ens = lh_ens[11:13,...].mean(axis=0)
    lh_ens = lh_ens[(lats<-7)&(lats>-47),:]
    lh_ens = lh_ens[:,(lons<155)&(lons>110)]
    lh_nino[j,...] = lh_ens

#lh_nina = np.ones((30,len(mlats),len(mlons)))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/'
#for j,ens in enumerate(lanina):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    lh_ens = np.squeeze(ncfile.variables['lh'][:,0,...])
#    lh_ens = lh_ens[11:13,...].mean(axis=0)
#    lh_ens = lh_ens[(lats<-7)&(lats>-47),:]
#    lh_ens = lh_ens[:,(lons<155)&(lons>110)]
#    lh_nina[j,...] = lh_ens
mp = Basemap(ax=axes[1],projection='cyl',llcrnrlon=110, llcrnrlat=-47, urcrnrlon=155, urcrnrlat=-7)
mp.drawmeridians([120,130,140,150],labels=[1,0,0,1],dashes=[5,700],fontsize=9)
xx,yy = mp(x,y)
dat = lh_nino[:11].mean(axis=0)-lh_nino[11:].mean(axis=0)
dat[mask==0] = 0
dat = np.ma.array(dat, mask=mask<50)
mesh = mp.pcolormesh(xx,yy, dat, cmap='PuOr',vmin=-15,vmax=15)
cont = mp.contour(xx,yy,dat,levels=levels,colors='k',linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
#_, p = stats.ttest_ind(lh_nino[:11], lh_nino[11:], axis=0, equal_var=False)
#mp.contourf(xx,yy,p<0.05,1,colors='none',hatches=[None,'xx'])
cbar = mp.colorbar(mesh,location='bottom',ticks=levels,pad=0.2)
cbar.set_label('Wm$^{-2}$')
mp.drawcoastlines()
plt.sca(axes[1])
plt.title('b) Q$_{E}$',loc='left')

#------------------------------------------------------

sm_nino = np.ones((30,len(mlats),len(mlons)))*np.nan
ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/'
for j,ens in enumerate(nino_few+nino_many):
    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
    sm_ens = np.squeeze(ncfile.variables['sm'][:,0,...])
    sm_ens = sm_ens[11:13,...].mean(axis=0)
    sm_ens = sm_ens[(lats<-7)&(lats>-47),:]
    sm_ens = sm_ens[:,(lons<155)&(lons>110)]
    sm_nino[j,...] = sm_ens

#sm_nina = np.ones((30,len(mlats),len(mlons)))*np.nan
#ensdir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/'
#for j,ens in enumerate(lanina):
#    ncfile = nc.MFDataset(ensdir+ens+'/'+ens+'a.pa2*')
#    sm_ens = np.squeeze(ncfile.variables['sm'][:,0,...])
#    sm_ens = sm_ens[11:13,...].mean(axis=0)
#    sm_ens = sm_ens[(lats<-7)&(lats>-47),:]
#    sm_ens = sm_ens[:,(lons<155)&(lons>110)]
#    sm_nina[j,...] = sm_ens
mp = Basemap(ax=axes[2],projection='cyl',llcrnrlon=110, llcrnrlat=-47, urcrnrlon=155, urcrnrlat=-7)
mp.drawmeridians([120,130,140,150],labels=[1,0,0,1],dashes=[5,700],fontsize=9)
xx,yy = mp(x,y)
dat = sm_nino[:11].mean(axis=0)-sm_nino[11:].mean(axis=0)
dat[mask==0] = 0
dat = np.ma.array(dat, mask=mask<50)
mesh = mp.pcolormesh(xx,yy,dat, cmap='BrBG',vmin=-0.5,vmax=0.5)
levels = np.arange(-0.5,0.6,0.2)
dat = np.ma.array(dat, mask=mask<50)
cont = mp.contour(xx,yy,dat,levels=levels,colors='k',linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
#_, p = stats.ttest_ind(np.ma.array(sm_nino[:11], mask=mask<50), np.ma.array(sm_nino[11:], mask=mask<50), axis=0, equal_var=False)
#mp.contourf(xx,yy,p<0.05,1,colors='none',hatches=[None,'xx'])
cbar = mp.colorbar(mesh,location='bottom',pad=0.2)
cbar.set_label('kgm$^{-2}$')
cbar.set_ticks(levels)
mp.drawcoastlines()
#mp.drawparallels([-10,-20,-30,-40],lables=[1,0,0,1],linewidth=0.03,fontsize=10)
plt.sca(axes[2])
plt.title('c) Soil Moisture',loc='left')
#plt.show()
plt.savefig('nino_soil.png',dpi=250,format='png')

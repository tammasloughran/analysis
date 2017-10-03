# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 16:17:33 2017

@author: tammas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import pandas as pd
import glob
import scipy.stats as stats
import pdb
import sys

print "Loading control data"
# Define directory and ensembles
data_directory = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
pacnino_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/pac_nino/va*')
pacnina_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/pac_nina/va*')
indpiod_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/ind_piod/va*')
indniod_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/ind_niod/va*')
indpac_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/indpac_nino/va*')
control = []
for year in xrange(1969,2000):
    for month in xrange(1,13):
        control += [data_directory+'vamrb/vamrba.pa'+str(year)+'-'+"%02d"%(month,)+'.nc']

# load climatology from control
controlnc = nc.MFDataset(control)
times = controlnc.variables['t'][:]
dates = nc.num2date(times,units=controlnc.variables['t'].units)
date_range = pd.period_range(dates[0],dates[-1], freq='M')
summer = (date_range.month==12)|(date_range.month==1)|(date_range.month==2)
# Pressure
mslp_hold = np.squeeze(controlnc.variables['p'][:])
mslp_hold = mslp_hold[summer]
mslp_ctrl = np.zeros((30,)+mslp_hold.shape[1:])
for yr in xrange(0,30):
    mslp_ctrl[yr,...] = mslp_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
mslp_summer_clim = np.mean(mslp_ctrl, axis=0)
# Temperature
temp_hold = np.squeeze(controlnc.variables['temp'][:])
temp_hold = temp_hold[summer]
temp_ctrl = np.zeros((30,)+temp_hold.shape[1:])
for yr in xrange(0,30):
    temp_ctrl[yr,...] = temp_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
temp_summer_clim = np.mean(temp_ctrl, axis=0)
# Rain
rain_hold = np.squeeze(controlnc.variables['precip'][:])*86400 # *86400 converts from kg m-2 s-1 to mm day-1
rain_hold = rain_hold[summer]
rain_ctrl = np.zeros((30,)+rain_hold.shape[1:])
for yr in xrange(0,30):
    rain_ctrl[yr,...] = rain_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
rain_summer_clim = np.mean(rain_ctrl, axis=0)
# wind
#u_ctrl = controlnc.variables['u'][:]
#u_ctrl = u_ctrl[summer]
#u_summer_clim = np.squeeze(np.mean(u_ctrl, axis=0))
#v_ctrl = controlnc.variables['u'][:]
#v_ctrl = v_ctrl[summer]
#v_summer_clim = np.squeeze(np.mean(v_ctrl, axis=0))
# 200hPa geopotential
pl = 200
p_1 = controlnc.variables['p_1'][:]
z = np.where(p_1==pl)[0][0]
hgt_hold = np.squeeze(controlnc.variables['ht_1'][:,z,...])
hgt_hold = hgt_hold[summer]
hgt_ctrl = np.zeros((30,)+hgt_hold.shape[1:])
for yr in xrange(0,30):
    hgt_ctrl[yr,...] = hgt_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
hgt_summer_clim = np.mean(hgt_ctrl, axis=0)

# Load latitude and longitude
lats = controlnc.variables['latitude'][:]
lats1 = controlnc.variables['latitude_1'][:]
lons = controlnc.variables['longitude'][:]


# Function to loop and mean for all ensembles
def get_data(experiment_ens):
    mslp = np.zeros((30,)+mslp_summer_clim.shape)
    temp = np.zeros((30,)+temp_summer_clim.shape)
    rain = np.zeros((30,)+rain_summer_clim.shape)
    #u = np.zeros((30,)+u_summer_clim.shape)
    #v = np.zeros((30,)+v_summer_clim.shape)
    hgt = np.zeros((30,)+hgt_summer_clim.shape)
    for n,ens in enumerate(experiment_ens):
        print ens
        summer_files = glob.glob(ens+'/*.pa2000-12.nc') # December
        summer_files += glob.glob(ens+'/*.pa2001-01.nc') # Janurary
        summer_files += glob.glob(ens+'/*.pa2001-02.nc') # February
        summernc = nc.MFDataset(summer_files)
        mslp[n,...] = np.squeeze(np.mean(summernc.variables['p'][...], axis=0))
        temp[n,...] = np.squeeze(np.mean(summernc.variables['temp'][:], axis=0))
        rain[n,...] = np.squeeze(np.mean(summernc.variables['precip'][:], axis=0))*86400
        #u_ensmean = np.squeeze(np.mean(summernc.variables['u'][:], axis=0))
        #v_ensmean = np.squeeze(np.mean(summernc.variables['v'][:], axis=0))
        hgt[n,...] = np.squeeze(np.mean(summernc.variables['ht_1'][:,z,...], axis=0))
        summernc.close()
    return mslp, temp, rain, hgt


# Define plotting functions
def plot_map(data,sig,lats,lons,units,rng,fname):
    plt.clf()
    mp = Basemap(projection='robin', lon_0=180)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    if units!='$mm/day$':
        colours = 'bwr'
        plot_cont = True
    else: 
        colours = 'BrBG'
        plot_cont = True
    clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
    if plot_cont:
        step = 2.*rng/10.
        levels = np.arange(-rng,rng+step,step)
        mp.contour(xx,yy,data.squeeze(),colors='k',levels=levels,linewidths=0.5)
    mp.contourf(xx,yy,sig,1,colors='none',hatches=[None,'xx'])
    cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
    cb.set_label(units)
    mp.drawcoastlines()
    #pll = mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0.01)
    #for line in pll.iteritems(): line[1][0].pop(0)
    #mer = mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0.01)
    #for line in mer.iteritems(): line[1][0].pop(0)
    plt.savefig(fname,format='png',dpi=250)


def plot_aus(data,sig,lats,lons,units,rng,fname):
    plt.clf()
    mp = Basemap(projection='mill',
             llcrnrlon=110.,llcrnrlat=-48.,
             urcrnrlon=157.,urcrnrlat=-5.)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    if units!='$mm/day$':
        colours = 'bwr'
        plot_cont = True
    else: 
        colours = 'BrBG'
        plot_cont = True
    clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
    if plot_cont:
        step = 2.*rng/10.
        levels = np.arange(-rng,rng+step,step)
        mp.contour(xx,yy,data.squeeze(),colors='k',levels=levels,linewidths=0.5)
    mp.contourf(xx,yy,sig,1,colors='none',hatches=[None,'xx'])
    cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
    cb.set_label(units)
    mp.drawcoastlines()
    plt.savefig(fname,format='png',dpi=250)

def plot_stream(u,v,lats,lons,units,rng,fname):
    plt.clf()
    mp = Basemap(projection='cyl', lon_0=180.)
    velo = np.sqrt(np.square(u)+np.square(v))
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    mp.pcolormesh(xx,yy,velo,cmap='viridis',vmin=0,vmax=rng)
    mp.streamplot(xx,yy,u,v, color='k')
    cb = mp.colorbar(location='bottom',pad=0.3)
    mp.drawcoastlines()
    mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0.01)
    mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0.01)
    cb.set_label(units)
    plt.savefig(fname,format='png',dpi=250)

def signif(data1,data2):
    _, p = stats.ttest_ind(data1, data2, axis=0, equal_var=False, nan_policy='propagate')
    return np.squeeze(p<0.05)

print "Meaning Pacific El nino"
mslp, temp, rain, hgt = get_data(pacnino_ens)
sig = signif(mslp, mslp_ctrl)
plot_map(mslp.mean(axis=0)-mslp_summer_clim,sig,lats,lons,'Pa',pl,'mslp_pacnino_ensmean.png')
sig = signif(temp, temp_ctrl)
plot_map(temp.mean(axis=0)-temp_summer_clim,sig,lats,lons,'${}^{\circ}C$',3,'temp_pacnino_ensmean.png')
sig = signif(rain, rain_ctrl)
plot_aus(rain.mean(axis=0)-rain_summer_clim,sig,lats,lons,'$mm/day$', 7,'rain_pacnino_ensmean.png')
#plot_stream(u.mean(axis=0),v.mean(axis=0),sig,lats1,lons,'$ms^{-1}$',12,'winds_pacnino_ensmean.png')
sig = signif(hgt, hgt_ctrl)
plot_map(hgt.mean(axis=0)-hgt_summer_clim,sig,lats1,lons,'m', 80,'hgt_pacnino_ensmean.png')

print "Meaning Pacific La nina"
mslp, temp, rain, hgt = get_data(pacnina_ens)
sig = signif(mslp, mslp_ctrl)
plot_map(mslp.mean(axis=0)-mslp_summer_clim,sig,lats,lons,'Pa',pl,'mslp_pacnina_ensmean.png')
sig = signif(temp, temp_ctrl)
plot_map(temp.mean(axis=0)-temp_summer_clim,sig,lats,lons,'${}^{\circ}C$',3,'temp_pacnina_ensmean.png')
sig = signif(rain, rain_ctrl)
plot_map(rain.mean(axis=0)-rain_summer_clim,sig,lats,lons,'$mm/day$', 7,'rain_pacnina_ensmean.png')
#plot_stream(u.mean(axis=0),v.mean(axis=0),sig,lats1,lons,'$ms^{-1}$',12,'winds_pacnina_ensmean.png')
sig = signif(hgt, hgt_ctrl)
plot_map(hgt.mean(axis=0)-hgt_summer_clim,sig,lats1,lons,'m', 80,'hgt_pacnina_ensmean.png')

print "Meaning Indian PIOD"
mslp, temp, rain, hgt = get_data(indpiod_ens)
sig = signif(mslp, mslp_ctrl)
plot_map(mslp.mean(axis=0)-mslp_summer_clim,sig,lats,lons,'Pa',pl,'mslp_indpiod_ensmean.png')
sig = signif(temp, temp_ctrl)
plot_map(temp.mean(axis=0)-temp_summer_clim,sig,lats,lons,'${}^{\circ}C$',3,'temp_indpiod_ensmean.png')
sig = signif(rain, rain_ctrl)
plot_map(rain.mean(axis=0)-rain_summer_clim,sig,lats,lons,'$mm/day$', 7,'rain_indpiod_ensmean.png')
#plot_stream(u_piodmean,v_piodmean,sig,lats1,lons,'$ms^{-1}$',12,'winds_indpiod_ensmean.png')
sig = signif(hgt, hgt_ctrl)
plot_map(hgt.mean(axis=0)-hgt_summer_clim,sig,lats1,lons,'m', 80,'hgt_indpiod_ensmean.png')

print "Meaning Indian NIOD"
mslp, temp, rain, hgt = get_data(indniod_ens)
sig = signif(mslp, mslp_ctrl)
plot_map(mslp.mean(axis=0)-mslp_summer_clim,sig,lats,lons,'Pa',pl,'mslp_indniod_ensmean.png')
sig = signif(temp, temp_ctrl)
plot_map(temp.mean(axis=0)-temp_summer_clim,sig,lats,lons,'${}^{\circ}C$',3,'temp_indniod_ensmean.png')
sig = signif(rain, rain_ctrl)
plot_map(rain.mean(axis=0)-rain_summer_clim,sig,lats,lons,'$mm/day$', 7,'rain_indniod_ensmean.png')
#plot_stream(u.mean(axis=0),v.mean(axis=0),sig,lats1,lons,'$ms^{-1}$',12,'winds_indniod_ensmean.png')
sig = signif(hgt, hgt_ctrl)
plot_map(hgt.mean(axis=0)-hgt_summer_clim,sig,lats1,lons,'m', 80,'hgt_indniod_ensmean.png')

print "Meaning Indpac Nino"
mslp, temp, rain, hgt = get_data(indpac_ens)
sig = signif(mslp, mslp_ctrl)
plot_map(mslp.mean(axis=0)-mslp_summer_clim,sig,lats,lons,'Pa',pl,'mslp_indpac_ensmean.png')
sig = signif(temp, temp_ctrl)
plot_map(temp.mean(axis=0)-temp_summer_clim,sig,lats,lons,'${}^{\circ}C$',3,'temp_indpac_ensmean.png')
sig = signif(rain, rain_ctrl)
plot_map(rain.mean(axis=0)-rain_summer_clim,sig,lats,lons,'$mm/day$', 7,'rain_indpac_ensmean.png')
#plot_stream(u.mean(axis=0),v.mean(axis=0),sig,lats1,lons,'$ms^{-1}$',12,'winds_indpac_ensmean.png')
sig = signif(hgt, hgt_ctrl)
plot_map(hgt.mean(axis=0)-hgt_summer_clim,sig,lats1,lons,'m', 80,'hgt_indpac_ensmean.png')
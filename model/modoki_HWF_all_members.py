# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:09:19 2017

@author: Tammas Loughran

Plot all ensemble members for EP En Nino and Modoki Nino
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap


# Define directories and ensembles
hwdir = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
hwdir2 = '/srv/ccrc/data48/z5032520/ehfheatwaves/'
modoki = ['vaqoc','vaqog','vaqok','vaqoo','vaqos','vaqow','vaqpa','vaqod',
          'vaqoh','vaqol','vaqop','vaqot','vaqox','vaqpb','vaqoa','vaqoe',
          'vaqoi','vaqom','vaqoq','vaqou','vaqoy','vaqpc','vaqob','vaqof',
          'vaqoj','vaqon','vaqor','vaqov','vaqoz','vaqpd']
elnino = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg',
          'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq',
          'vaoqr','vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy',
          'vaoqz','vaqgl','vaqgm','vaqgn','vaoqh','vaoqn']
control = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh',
           'vaowf','vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo',
           'vaowp','vaowq','vaowr','vaows','vaowt','vaowu','vaowv','vaoww',
           'vaowx','vaowy','vaowz','vaqgi','vaqgj','vaqgk']
lanina = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl',
          'vamrm','vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt',
          'vamru','vamrv','vamrw','vamrx','vamry','vamrz','vaqga','vaqgb',
          'vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']

def load_hwf(filename):
    ncfile = nc.Dataset(filename, 'r')
    return ncfile.variables['HWF_EHF'][0,...]


# Load data from test file
filename = hwdir+'vamrd/EHF_heatwaves_ACCESS1.3_vamrd_yearly_summer.nc'
test_hwf = load_hwf(filename)
testnc = nc.Dataset(filename)
lats = testnc.variables['lat'][:]
lons = testnc.variables['lon'][:]
testnc.close()
x, y = np.meshgrid(lons, lats)
del lats, lons

# Control ensemble mean
control_hwf = np.ma.ones((30,)+test_hwf.shape)*np.nan
for i, ens in enumerate(control):
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    control_hwf[i,...] = load_hwf(filename)
control_mean_hwf = np.ma.mean(control_hwf, axis=0)

# Load data and express as anomalies from control mean.
elnino_hwf = np.ma.ones((30,)+test_hwf.shape)*np.nan
for i, ens in enumerate(elnino):
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    elnino_hwf[i,...] = load_hwf(filename) - control_mean_hwf

modoki_hwf = np.ma.ones((30,)+test_hwf.shape)*np.nan
for i, ens in enumerate(modoki):
    filename = hwdir2+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    modoki_hwf[i,...] = load_hwf(filename) - control_mean_hwf

lanina_hwf = np.ma.ones((30,)+test_hwf.shape)*np.nan
for i, ens in enumerate(lanina):
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    lanina_hwf[i,...] = load_hwf(filename) - control_mean_hwf

# Plot all the control ensembles
f, axes = plt.subplots(nrows=6, ncols=5, figsize=(10,9))
for i, ax in enumerate(axes.flat):
    lvls = np.arange(0,30,5)
    mp = Basemap(ax=ax, projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    xx, yy = mp(x,y)
    shade = mp.pcolormesh(xx,yy,control_hwf[i,...],
                          cmap='viridis',vmin=0,vmax=30)
    cont = mp.contour(xx,yy,control_hwf[i,...],levels=lvls,colors='k')
    for c in cont.collections:
        c.set_linewidth(0.5)
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    mp.drawcoastlines()
    ax.set_title(str(i+1)+')', loc='left', fontsize=10)
cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='vertical')
cb = plt.colorbar(shade, cax=cbar_ax, orientation='vertical', ticks=lvls)
ylbl = cb.ax.set_ylabel('days')
ylbl.set_rotation(0)
plt.show()

# Plot all the La Nina ensembles
f, axes = plt.subplots(nrows=6, ncols=5, figsize=(10,9))
for i, ax in enumerate(axes.flat):
    lvls = np.arange(-25,26,5)
    mp = Basemap(ax=ax, projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    xx, yy = mp(x,y)
    shade = mp.pcolormesh(xx,yy,lanina_hwf[i,...],
                          cmap='bwr',vmin=-25,vmax=25)
    cont = mp.contour(xx,yy,lanina_hwf[i,...],levels=lvls,colors='k')
    for c in cont.collections:
        c.set_linewidth(0.5)
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    mp.drawcoastlines()
    ax.set_title(str(i+1)+')', loc='left', fontsize=10)
cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='vertical')
cb = plt.colorbar(shade, cax=cbar_ax, orientation='vertical', ticks=lvls)
ylbl = cb.ax.set_ylabel('days')
ylbl.set_rotation(0)
plt.show()

# Plot all the EP El Nino ensembles
f, axes = plt.subplots(nrows=6, ncols=5, figsize=(10,9))
for i, ax in enumerate(axes.flat):
    lvls = np.arange(-25,26,5)
    mp = Basemap(ax=ax, projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    xx, yy = mp(x,y)
    shade = mp.pcolormesh(xx,yy,elnino_hwf[i,...],
                          cmap='bwr',vmin=-25,vmax=25)
    cont = mp.contour(xx,yy,elnino_hwf[i,...],levels=lvls,colors='k',linewidth=0.5)
    for c in cont.collections:
        c.set_linewidth(0.5)
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    mp.drawcoastlines()
    ax.set_title(str(i+1)+')', loc='left', fontsize=10)
cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='vertical')
cb = plt.colorbar(shade, cax=cbar_ax, orientation='vertical', ticks=lvls)
ylbl = cb.ax.set_ylabel('days')
ylbl.set_rotation(0)
plt.show()

# Plot all the Modoki ensembles
f, axes = plt.subplots(nrows=6, ncols=5, figsize=(10,9))
for i, ax in enumerate(axes.flat):
    lvls = np.arange(-25,26,5)
    mp = Basemap(ax=ax, projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    xx, yy = mp(x,y)
    shade = mp.pcolormesh(xx,yy,modoki_hwf[i,...],
                          cmap='bwr',vmin=-25,vmax=25)
    cont = mp.contour(xx,yy,modoki_hwf[i,...],levels=lvls,colors='k',linewidth=0.5)
    for c in cont.collections:
        c.set_linewidth(0.5)
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    mp.drawcoastlines()
    ax.set_title(str(i+1)+')', loc='left', fontsize=10)
cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='vertical')
cb = plt.colorbar(shade, cax=cbar_ax, orientation='vertical', ticks=lvls)
ylbl = cb.ax.set_ylabel('days')
ylbl.set_rotation(0)
plt.show()
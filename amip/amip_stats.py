# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:55:12 2017

@author: Tammas Loughran
This scripts plots some basic multi-model and individual-model
heatwave statistics from the AMIP experiment of CMIP5. 
"""
import os
import sys
import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap

# User. Manually select the aspect and enso phase
aspect = 'HWM'
phase = 'lanina'

# Define some variables
cwd = os.getcwd()
amipdir = '/srv/ccrc/data48/z5032520/amip/'
models = os.listdir(amipdir)
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CMCC-CM','CNRM-CM5',
          'CSIRO-Mk3-6-0','GFDL-CM3','HadGEM2-A','IPSL-CM5A-LR','MIROC5',
          'MPI-ESM-MR','MRI-CGCM3','NorESM1-M','inmcm4']
#models = ['IPSL-CM5A-LR','MIROC5','MPI-ESM-MR','MRI-CGCM3','NorESM1-M','inmcm4']
elnino_years = [1982, 1987, 1991, 1997, 2002]
lanina_years = [1984, 1988, 1998, 1999, 2007]
if phase=='elnino':
    phase_years = elnino_years
elif phase=='lanina':
    phase_years = lanina_years

parallels = np.arange(-40., -9., 10.)
meridians = np.arange(120., 160., 10.,)


def plot_map(flons, flats, data, title='unknown'):
    """Plot a map of heatwave data
    """
    if 'HWF' in title:
        levels = np.arange(0.,21.,2.)
    elif 'HWD' in title:
        levels = np.arange(0.,11.,2.)
    elif 'HWN' in title:
        levels = np.arange(0.,5.5,1.)
    elif 'HWA' in title:
        levels = np.arange(0.,51.,5.)
    elif 'HWM' in title:
        levels = np.arange(0.,31.,5.)
    #rng = np.arange(data.min(),data.max(),2)
    m = Basemap(projection='mill',
            llcrnrlon=110.,llcrnrlat=-45.,
            urcrnrlon=157.,urcrnrlat=-10.,
            resolution='l')
    lns, lts = np.meshgrid(flons, flats)
    x,y = m(lns,lts)
    shade = m.pcolormesh(x,y,data,cmap='YlOrRd',vmin=0,vmax=levels[-1])
    #cont = m.contour(x,y,data,linewidths=0.4,colors='k',levels=levels)
    #for c in cont.collections:
    #    if c.get_linestyle() == [(None, None)]:
    #       continue
    #    else:
    #        c.set_dashes([(0, (2.0, 2.0))])
    m.drawcoastlines()
    m.drawmeridians(meridians, linewidth=0, labels=[False,False,False,True])
    m.drawparallels(parallels, linewidth=0, labels=[True,False,False,False])
    plt.title(title)
    m.colorbar(shade)
    plt.savefig(title+'.eps', format='eps')
    plt.close()

def plot_map_div(flons, flats, data, ax=None, title='unknown'):
    """Plot a map of heatwave data
    """
    if 'HWF' in title:
        levels = np.arange(-20.,21.,2.)
    elif 'HWD' in title:
        levels = np.arange(-8.,9.,1.)
    elif 'HWN' in title:
        levels = np.arange(-3.,3.5,0.5)
    elif 'HWA' in title:
        levels = np.arange(-30.,31.,5.)
    elif 'HWM' in title:
        levels = np.arange(-16.,17.,2.)
    m = Basemap(ax=ax, projection='mill',
            llcrnrlon=110.,llcrnrlat=-45.,
            urcrnrlon=157.,urcrnrlat=-10.,
            resolution='l')
    lns, lts = np.meshgrid(flons, flats)
    x,y = m(lns,lts)
    shade = m.pcolormesh(x,y,data,cmap='bwr',vmin=levels[0],vmax=levels[-1])
    cont = m.contour(x,y,data,colors='k',linewidths=0.4,levels=levels[::2])
    for c in cont.collections:
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    m.drawcoastlines()
    ax.annotate(title,(0.02,0.02),xycoords='axes fraction')
    return shade

def locate_in_list(values,alist):
    """Locates the indices of 'alist' where the first instance of values
    in 'values' occurrs.
    
    Input
    values - a list of values to locate
    alist - the list to search
    
    Returns
    index - a list of indices
    """
    index = []
    for i in values:
        index.append(np.where(alist==i)[0][0])
    return index


# Load AWAP Data
awapnc = nc.Dataset('/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/EHF_heatwaves_AWAP_bp1979-2008_yearly_summer.nc')
var = awapnc.variables[aspect+'_EHF'][...]
lons = awapnc.variables['lon'][:]
lats = awapnc.variables['lat'][:]
times = awapnc.variables['time'][:]
phase_i = locate_in_list(phase_years, times)
clim = (times>=1979)&(times<=2008)
aspect_phase = var[phase_i,...].mean(axis=0)-var[clim].mean(axis=0)

# Construct figure
f, axes = plt.subplots(nrows=4, ncols=4, figsize=(10,9), sharex=True, sharey=True)
f.subplots_adjust(hspace=0, wspace=0)
# Plot AWAP
mappable = plot_map_div(lons, lats, aspect_phase, ax=axes[0][0], title=aspect+'_AWAP')


for model,ax in zip(models,[x for x in axes.flat][1:]):
    print model
    
    # Load data
    hwfiles = glob.glob(amipdir+model+'/EHF*yearly*')
    dataset = nc.Dataset(hwfiles[0])
    shape = dataset.variables[aspect+'_EHF'].shape
    ifile = hwfiles[0]
    hwnc = nc.Dataset(ifile)
    var = hwnc.variables[aspect+'_EHF'][...]
    lons = hwnc.variables['lon'][:]
    lons = lons - (lons[1] - lons[0])/2.
    lats = hwnc.variables['lat'][:]
    lats = lats - (lats[1] - lats[0])/2.
    time = hwnc.variables['time'][:]
    clim = (time>=1979)&(time<=2008)
    phase_i = locate_in_list(phase_years, time)
    
    # Construct the compsites (expressed as anomaly from climatology)
    aspect_phase = var[phase_i,...].mean(axis=0)-var[clim].mean(axis=0)

    # Plot the elnino composites 
    plot_map_div(lons,lats,aspect_phase,ax=ax,title=aspect+'_'+model)


# Plot the colorbar
if aspect=='HWF':
    levels = np.arange(-20.,21.,2.)
elif aspect=='HWD':
    levels = np.arange(-8.,9.,1.)
elif aspect=='HWN':
    levels = np.arange(-3.,3.5,0.5)
elif aspect=='HWA':
    levels = np.arange(-30.,31.,5.)
elif aspect=='HWM':
    levels = np.arange(-16.,17.,2.)
cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='horizontal',pad=0.02)
cb = plt.colorbar(mappable, cax=cbar_ax, orientation='horizontal', ticks=levels[::2])

# Save figure
plt.savefig('AMIP_'+aspect+'_'+phase+'.eps',format='eps')
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:55:12 2017

@author: Tammas Loughran
This scripts plots some basic multi-model and individual-model
heatwave statistics from the AMIP experiment of CMIP5. 
"""
import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


# Define some variables
cwd = os.getcwd()
amipdir = '/srv/ccrc/data48/z5032520/amip/'
models = os.listdir(amipdir)
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CMCC-CM','CNRM-CM5',
          'CSIRO-Mk3-6-0','GFDL-CM3','HadGEM2-A','IPSL-CM5A-LR','MIROC5',
          'MPI-ESM-MR','MRI-CGCM3','NorESM1-M','inmcm4']
models = ['IPSL-CM5A-LR','MIROC5','MPI-ESM-MR','MRI-CGCM3','NorESM1-M','inmcm4']
elnino_years = [1982, 1987, 1991, 1997, 2002]
lanina_years = [1984, 1988, 1998, 1999, 2007]
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

def plot_map_div(flons, flats, data, title='unknown'):
    """Plot a map of heatwave data
    """
    if 'HWF' in title:
        levels = np.arange(-20.,21.,2.)
    elif 'HWD' in title:
        levels = np.arange(-10.,11.,2.)
    elif 'HWN' in title:
        levels = np.arange(-5.,5.5,1.)
    elif 'HWA' in title:
        levels = np.arange(-50.,51.,5.)
    elif 'HWM' in title:
        levels = np.arange(-30.,31.,5.)
    #vmax = round(max(data.max(),abs(data.min())))
    m = Basemap(projection='mill',
            llcrnrlon=110.,llcrnrlat=-45.,
            urcrnrlon=157.,urcrnrlat=-10.,
            resolution='l')
    lns, lts = np.meshgrid(flons, flats)
    x,y = m(lns,lts)
    shade = m.pcolormesh(x,y,data,cmap='bwr',vmin=levels[0],vmax=levels[-1])
    #cont = m.contour(x,y,data,colors='k',linewidths=0.4,levels=levels)
    m.drawcoastlines()
    m.drawmeridians(meridians, linewidth=0, labels=[False,False,False,True])
    m.drawparallels(parallels, linewidth=0, labels=[True,False,False,False])
    plt.title(title)
    m.colorbar(shade)
    plt.savefig(title+'.eps', format='eps')
    plt.close()

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


for model in models:
    print model
    # Make model directory
    os.chdir(cwd)
    try: os.mkdir(model)
    except: pass
    os.chdir(model)
    
    # Load data
    hwfiles = glob.glob(amipdir+model+'/EHF*yearly*')
    dataset = nc.Dataset(hwfiles[0])
    shape = dataset.variables['HWF_EHF'].shape
    ens = 1
    hwf_ens_mean = np.zeros(shape[1:])
    hwd_ens_mean = np.zeros(shape[1:])
    hwn_ens_mean = np.zeros(shape[1:])
    hwa_ens_mean = np.zeros(shape[1:])
    hwm_ens_mean = np.zeros(shape[1:])
    hwf_ens_elnino = np.zeros(shape[1:])
    hwd_ens_elnino = np.zeros(shape[1:])
    hwn_ens_elnino = np.zeros(shape[1:])
    hwa_ens_elnino = np.zeros(shape[1:])
    hwm_ens_elnino = np.zeros(shape[1:])
    hwf_ens_lanina = np.zeros(shape[1:])
    hwd_ens_lanina = np.zeros(shape[1:])
    hwn_ens_lanina = np.zeros(shape[1:])
    hwa_ens_lanina = np.zeros(shape[1:])
    hwm_ens_lanina = np.zeros(shape[1:])
    for ifile in hwfiles:
        print ens
        hwnc = nc.Dataset(ifile)
        hwf = hwnc.variables['HWF_EHF'][...]
        hwd = hwnc.variables['HWD_EHF'][...]
        hwn = hwnc.variables['HWN_EHF'][...]
        hwa = hwnc.variables['HWA_EHF'][...]
        hwm = hwnc.variables['HWM_EHF'][...]
        lons = hwnc.variables['lon'][:]
        lons = lons - (lons[1] - lons[0])/2.
        lats = hwnc.variables['lat'][:]
        lats = lats - (lats[1] - lats[0])/2.
        time = hwnc.variables['time'][:]
        elnino_i = locate_in_list(elnino_years, time)
        lanina_i = locate_in_list(lanina_years, time)    
        
        # Plot the climatology.
        plot_map(lons,lats,hwf.mean(axis=0),title='HWF_mean_'+model+'_'+str(ens))
        plot_map(lons,lats,hwd.mean(axis=0),title='HWD_mean_'+model+'_'+str(ens))
        plot_map(lons,lats,hwn.mean(axis=0),title='HWN_mean_'+model+'_'+str(ens))
        plot_map(lons,lats,hwa.mean(axis=0),title='HWA_mean_'+model+'_'+str(ens))
        plot_map(lons,lats,hwm.mean(axis=0),title='HWM_mean_'+model+'_'+str(ens))
        
        # Plot the standard deviation
        plot_map(lons,lats,np.std(hwf,axis=0),title='HWF_std_'+model+'_'+str(ens))
        plot_map(lons,lats,np.std(hwd,axis=0),title='HWD_std_'+model+'_'+str(ens))
        plot_map(lons,lats,np.std(hwn,axis=0),title='HWN_std_'+model+'_'+str(ens))
        plot_map(lons,lats,np.std(hwa,axis=0),title='HWA_std_'+model+'_'+str(ens))
        plot_map(lons,lats,np.std(hwm,axis=0),title='HWM_std_'+model+'_'+str(ens))
        
        # Construct the compsites (expressed as anomaly from climatology)
        hwf_elnino = hwf[elnino_i,...].mean(axis=0)-hwf.mean(axis=0)
        hwd_elnino = hwd[elnino_i,...].mean(axis=0)-hwd.mean(axis=0)
        hwn_elnino = hwn[elnino_i,...].mean(axis=0)-hwn.mean(axis=0)
        hwa_elnino = hwa[elnino_i,...].mean(axis=0)-hwa.mean(axis=0)
        hwm_elnino = hwm[elnino_i,...].mean(axis=0)-hwm.mean(axis=0)
        hwf_lanina = hwf[lanina_i,...].mean(axis=0)-hwf.mean(axis=0)
        hwd_lanina = hwd[lanina_i,...].mean(axis=0)-hwd.mean(axis=0)
        hwn_lanina = hwn[lanina_i,...].mean(axis=0)-hwn.mean(axis=0)
        hwa_lanina = hwa[lanina_i,...].mean(axis=0)-hwa.mean(axis=0)
        hwm_lanina = hwm[lanina_i,...].mean(axis=0)-hwm.mean(axis=0)

        # Plot the elnino composites 
        plot_map_div(lons,lats,hwf_elnino,title='HWF_elnino_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwd_elnino,title='HWD_elnino_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwn_elnino,title='HWN_elnino_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwa_elnino,title='HWA_elnino_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwm_elnino,title='HWM_elnino_'+model+'_'+str(ens))
        
        # Plot the lanina composites
        plot_map_div(lons,lats,hwf_lanina,title='HWF_lanina_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwd_lanina,title='HWD_lanina_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwn_lanina,title='HWN_lanina_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwa_lanina,title='HWA_lanina_'+model+'_'+str(ens))
        plot_map_div(lons,lats,hwm_lanina,title='HWM_lanina_'+model+'_'+str(ens))
        
        # Aggregate Ensemble means      
        hwf_ens_mean += hwf.mean(axis=0)
        hwd_ens_mean += hwd.mean(axis=0)
        hwn_ens_mean += hwn.mean(axis=0)
        hwm_ens_mean += hwa.mean(axis=0)
        hwa_ens_mean += hwm.mean(axis=0)
        hwf_ens_elnino += hwf_elnino
        hwd_ens_elnino += hwd_elnino
        hwn_ens_elnino += hwn_elnino
        hwa_ens_elnino += hwa_elnino
        hwm_ens_elnino += hwm_elnino
        hwf_ens_lanina += hwf_lanina
        hwd_ens_lanina += hwd_lanina
        hwn_ens_lanina += hwn_lanina
        hwa_ens_lanina += hwa_lanina
        hwm_ens_lanina += hwm_lanina
        
        ens += 1
    
    # Divide the ensemble mean
    hwf_ens_mean = np.ma.array(hwf_ens_mean/float(ens-1),mask=hwf.mask[0])
    hwd_ens_mean = np.ma.array(hwd_ens_mean/float(ens-1),mask=hwf.mask[0])
    hwn_ens_mean = np.ma.array(hwn_ens_mean/float(ens-1),mask=hwf.mask[0])
    hwa_ens_mean = np.ma.array(hwa_ens_mean/float(ens-1),mask=hwf.mask[0])
    hwm_ens_mean = np.ma.array(hwm_ens_mean/float(ens-1),mask=hwf.mask[0])
    hwf_ens_elnino = np.ma.array(hwf_ens_elnino/float(ens-1),mask=hwf.mask[0])
    hwd_ens_elnino = np.ma.array(hwd_ens_elnino/float(ens-1),mask=hwf.mask[0])
    hwn_ens_elnino = np.ma.array(hwn_ens_elnino/float(ens-1),mask=hwf.mask[0])
    hwa_ens_elnino = np.ma.array(hwa_ens_elnino/float(ens-1),mask=hwf.mask[0])
    hwm_ens_elnino = np.ma.array(hwm_ens_elnino/float(ens-1),mask=hwf.mask[0])
    hwf_ens_lanina = np.ma.array(hwf_ens_lanina/float(ens-1),mask=hwf.mask[0])
    hwd_ens_lanina = np.ma.array(hwd_ens_lanina/float(ens-1),mask=hwf.mask[0])
    hwn_ens_lanina = np.ma.array(hwn_ens_lanina/float(ens-1),mask=hwf.mask[0])
    hwa_ens_lanina = np.ma.array(hwa_ens_lanina/float(ens-1),mask=hwf.mask[0])
    hwm_ens_lanina = np.ma.array(hwm_ens_lanina/float(ens-1),mask=hwf.mask[0])
  
    # Plot ensemble mean
    plot_map(lons,lats,hwf_ens_mean,title='HWF_mean_'+model+'_ensmean')
    plot_map(lons,lats,hwd_ens_mean,title='HWD_mean_'+model+'_ensmean')
    plot_map(lons,lats,hwn_ens_mean,title='HWN_mean_'+model+'_ensmean')
    plot_map(lons,lats,hwa_ens_mean,title='HWA_mean_'+model+'_ensmean')
    plot_map(lons,lats,hwm_ens_mean,title='HWM_mean_'+model+'_ensmean')
    plot_map_div(lons,lats,hwf_ens_elnino,title='HWF_elnino_'+model+'_ensmean')
    plot_map_div(lons,lats,hwd_ens_elnino,title='HWD_elnino_'+model+'_ensmean')
    plot_map_div(lons,lats,hwn_ens_elnino,title='HWN_elnino_'+model+'_ensmean')
    plot_map_div(lons,lats,hwa_ens_elnino,title='HWA_elnino_'+model+'_ensmean')
    plot_map_div(lons,lats,hwm_ens_elnino,title='HWM_elnino_'+model+'_ensmean')
    plot_map_div(lons,lats,hwf_ens_lanina,title='HWF_lanina_'+model+'_ensmean')
    plot_map_div(lons,lats,hwd_ens_lanina,title='HWD_lanina_'+model+'_ensmean')
    plot_map_div(lons,lats,hwn_ens_lanina,title='HWN_lanina_'+model+'_ensmean')
    plot_map_div(lons,lats,hwa_ens_lanina,title='HWA_lanina_'+model+'_ensmean')
    plot_map_div(lons,lats,hwm_ens_lanina,title='HWM_lanina_'+model+'_ensmean')
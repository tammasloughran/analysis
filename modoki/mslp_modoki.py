# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 14:57:22 2017

@author: Tammas Loughran

mslp_modoki.py
Plots the MSLPA and SSTA associated with modoki for both obs and ACCESS s
imulations.

"""
import netCDF4 as nc
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import glob

# Ensemble codes
modoki = ['vaqoc','vaqog','vaqok','vaqoo','vaqos','vaqow','vaqpa','vaqod',
          'vaqoh','vaqol','vaqop','vaqot','vaqox','vaqpb','vaqoa','vaqoe',
          'vaqoi','vaqom','vaqoq','vaqou','vaqoy','vaqpc','vaqob','vaqof',
          'vaqoj','vaqon','vaqor','vaqov','vaqoz','vaqpd']
elnino = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg',
          'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq',
          'vaoqr','vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy',
          'vaoqz','vaqgl','vaqgm','vaqgn','vaoqh','vaoqn']

def load_ensemble_hw(filename, hwdefinition='EHF', get_latlon=False):
    """Load the Australian hw ensemble data and lats and lons.

    Arguments:
    filename -- the path and filename of the file to load.
    hwdefinition -- the heatwave definition the file contains.

    Returns:
    hwf -- frequency
    hwn -- number
    hwd -- duration
    hwa -- amplitude
    hwm -- magnitude
    hwt -- timing
    lats -- latitudes
    lons -- longitudes
    """
    ncfile = nc.Dataset(filename)

    lats = ncfile.variables['lat'][:]
    lons = ncfile.variables['lon'][:]

    hwf = ncfile.variables['HWF_'+hwdefinition][0]
    hwf = hwf[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    hwn = ncfile.variables['HWN_'+hwdefinition][0]
    hwn = hwn[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    hwd = ncfile.variables['HWD_'+hwdefinition][0]
    hwd = hwd[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    hwd.data[hwd.mask] = np.nan

    hwa = ncfile.variables['HWA_'+hwdefinition][0]
    hwa = hwa[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    hwa.data[hwa.mask] = np.nan

    hwm = ncfile.variables['HWM_'+hwdefinition][0]
    hwm = hwm[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    hwm.data[hwm.mask] = np.nan

    hwt = ncfile.variables['HWT_'+hwdefinition][0]
    hwt = hwt[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    hwt.data[hwt.mask] = np.nan

    lats = lats[(lats<-10.)&(lats>-44.)]
    lons = lons[(lons<156.)&(lons>112.)]

    if get_latlon==False:
        lats, lons = 0, 0

    return hwf, hwn, hwd, hwa, hwm, hwt, lats, lons





def load_mslp_ensemble(ensemble, exp):
    """load mean sea level pressure of an ensemble member
    """
    if exp=='elnino': 
        data = 'data46'
    elif exp=='modokielnino':
        data = 'data48'
    filename = '/srv/ccrc/'+data+'/z5032520/modelout/ACCESS/'+exp+'/'+ensemble+'/'+ensemble+'a.pe*.nc'
    ncfile = nc.MFDataset(filename,'r')
    #time = ncfile.variables['t'][:]
    #tunits = ncfile.variables['t'].units
    #start = nc.num2date(time[0],units=tunits)
    #end = nc.num2date(time[-1],units=tunits)
    date_range = pd.period_range(start='2000-01-01',end='2001-12-31',freq='D')
    nov1 = dt.datetime(2000,11,1)
    mar31 = dt.datetime(2001,3,31)
    nov2mar = (date_range>=nov1)&(date_range<=mar31)
    pressure = ncfile.variables['p_1'][nov2mar,...]
    return np.squeeze(pressure.mean(axis=0))


if __name__=='__main__':
    # Load the pressure climatology
    print "loading"
    ncs = nc.MFDataset('/srv/ccrc/data46/z5032520/modelout/ACCESS/vamrb/vamrba.pa19*')
    prc = np.squeeze(ncs.variables['p'][-12*30:])
    lons = ncs.variables['longitude'][:]
    lats = ncs.variables['latitude'][:]
    ncs.close()
    date_range = pd.period_range(start='1970-01-01',end='1999-12-31',freq='M')
    pclim = np.zeros((12,)+prc.shape[-2:])
    print "clim"
    for month in xrange(1,13,1):
        pclim[month-1,...] = prc[date_range.month==month].mean(axis=0)
    m = [1,2,3,4,5,6,7,8,9,10,11,12]
    ndjfm = (m==1)&(m==2)&(m==3)&(m==11)&(m==12)
    ndjfm_pclim = pclim[ndjfm]

    #calculate the modoki ensemble mean
    modoki_mean = np.zeros(ndjfm_pclim.shape)
    for ens in modoki:
        print ens
        ens_pr = load_mslp_ensemble(ens,'modokielnino')
        ens_pr = ens_pr - ndjfm_pclim
        modoki_mean += ens_pr
    modoki_mean /= len(modoki)

    # calculate the elnino ensemble mean 
    elnino_mean = np.zeros(ndjfm_pclim.shape)
    for ens in elnino:
        print ens
        ens_pr = load_mslp_ensemble(ens,'elnino')
        ens_pr = ens_pr - ndjfm_pclim
        elnino_mean += ens_pr
    elnino_mean /= len(elnino)
    
    f, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
    m = Basemap(ax=ax1, projection='mill',
                llcrnrlon=110.,llcrnrlat=-45.,
                urcrnrlon=157.,urcrnrlat=-10.,
                resolution='l',fix_aspect=True)
    lns, lts = np.meshgrid(lons, lats)
    x,y = m(lns,lts)
    xx,yy = m(lns-.5,lts-.5)
    cont = m.contourf(xx,yy,modoki_mean,cmap='seismic')
    plt.colorbar(cont)
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    

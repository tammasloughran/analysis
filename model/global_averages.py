# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 17:13:46 2017

@author: tammas
"""

import numpy as np
#import scipy.stats as stats
import netCDF4 as nc
import glob
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

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

print "Loading:"
# load climatology from control
controlnc = nc.MFDataset(control)
lats = controlnc.variables['latitude'][:]
lons = controlnc.variables['longitude'][:]
times = controlnc.variables['t'][:]
dates = nc.num2date(times,units=controlnc.variables['t'].units)
date_range = pd.period_range(dates[0],dates[-1], freq='M')
summer = (date_range.month==12)|(date_range.month==1)|(date_range.month==2)
print("control temp")
# temp
temp_hold = np.squeeze(controlnc.variables['temp'][:,...])
temp_hold = temp_hold[summer]
ctrl_temp = np.zeros((30,)+temp_hold.shape[1:])
for yr in xrange(0,30):
    ctrl_temp[yr,...] = temp_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
print("Control rh")
# rh
rh_hold = np.squeeze(controlnc.variables['rh'][:,...])
rh_hold = rh_hold[summer]
ctrl_rh = np.zeros((30,)+rh_hold.shape[1:])
for yr in xrange(0,30):
    ctrl_rh[yr,...] = rh_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
print("control cloud")
# cloud
field30_1_hold = np.squeeze(controlnc.variables['field30_1'][:,...])
temp_field30_1 = field30_1_hold[summer]
ctrl_field30_1 = np.zeros((30,)+field30_1_hold.shape[1:])
for yr in xrange(0,30):
    ctrl_field30_1[yr,...] = field30_1_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
#sm_summer_clim = np.mean(ctrl_sm, axis=0)



# Load the ensemble soil moisture data
def get_ens(experiment_ens,var):
    summer_files = glob.glob(experiment_ens[0]+'/*.pa2000-12.nc')
    summernc = nc.MFDataset(summer_files)
    dshape = summernc.variables[var].shape[2:]
    summernc.close()
    nens = len(experiment_ens)
    data = np.zeros((nens,)+dshape)
    for n,ens in enumerate(experiment_ens):
        print ens
        summer_files = glob.glob(ens+'/*.pa2000-12.nc') # December
        summer_files += glob.glob(ens+'/*.pa2001-01.nc') # Janurary
        summer_files += glob.glob(ens+'/*.pa2001-02.nc') # February
        summernc = nc.MFDataset(summer_files)
        data[n,...] = np.squeeze(np.mean(summernc.variables[var][:,0,...],axis=0))
        summernc.close()
    #sm = sm[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
    #data = np.ma.array(data,mask=sm==0)
    return data


xx,yy = np.meshgrid(lons,lats)
coslats = np.cos(np.deg2rad(yy))
def global_ave(data,wgts):
    out = np.zeros(data.shape[0])
    for i in xrange(data.shape[0]):
        out[i] = np.average(data[i,...],weights=wgts)
    return out

control_temp_ave = global_ave(ctrl_temp,coslats) 
control_rh_ave = global_ave(ctrl_rh,coslats)
control_field30_1_ave = global_ave(ctrl_field30_1,coslats)


pacnino_temp = get_ens(pacnino_ens,'temp')
pacnina_temp = get_ens(pacnina_ens,'temp')
indpiod_temp = get_ens(indpiod_ens,'temp')
indniod_temp = get_ens(indniod_ens,'temp')
indpac_temp = get_ens(indpac_ens,'temp')
pacnino_temp_ave = global_ave(pacnino_temp,coslats)
pacnina_temp_ave = global_ave(pacnina_temp,coslats)
indpiod_temp_ave = global_ave(indpiod_temp,coslats)
indniod_temp_ave = global_ave(indniod_temp,coslats)
indpac_temp_ave = global_ave(indpac_temp,coslats)

#pacnino_q_2 = get_ens(pacnino_ens,'q_2')
#pacnina_q_2 = get_ens(pacnina_ens,'q_2')
#indpiod_q_2 = get_ens(indpiod_ens,'q_2')
#indniod_q_2 = get_ens(indniod_ens,'q_2')
#indpac_q_2 = get_ens(indpac_ens,'q_2')
#pacnino_q_2_ave = global_ave(pacnino_q_2,coslats)
#pacnina_q_2_ave = global_ave(pacnina_q_2,coslats)
#indpiod_q_2_ave = global_ave(indpiod_q_2,coslats)
#indniod_q_2_ave = global_ave(indniod_q_2,coslats)
#indpac_q_2_ave = global_ave(indpac_q_2,coslats)

pacnino_rh = get_ens(pacnino_ens,'rh')
pacnina_rh = get_ens(pacnina_ens,'rh')
indpiod_rh = get_ens(indpiod_ens,'rh')
indniod_rh = get_ens(indniod_ens,'rh')
indpac_rh = get_ens(indpac_ens,'rh')
pacnino_rh_ave = global_ave(pacnino_rh,coslats)
pacnina_rh_ave = global_ave(pacnina_rh,coslats)
indpiod_rh_ave = global_ave(indpiod_rh,coslats)
indniod_rh_ave = global_ave(indniod_rh,coslats)
indpac_rh_ave = global_ave(indpac_rh,coslats)

pacnino_field30_1 = get_ens(pacnino_ens,'field30_1')
pacnina_field30_1 = get_ens(pacnina_ens,'field30_1')
indpiod_field30_1 = get_ens(indpiod_ens,'field30_1')
indniod_field30_1 = get_ens(indniod_ens,'field30_1')
indpac_field30_1 = get_ens(indpac_ens,'field30_1')
pacnino_field30_1_ave = global_ave(pacnino_field30_1,coslats)
pacnina_field30_1_ave = global_ave(pacnina_field30_1,coslats)
indpiod_field30_1_ave = global_ave(indpiod_field30_1,coslats)
indniod_field30_1_ave = global_ave(indniod_field30_1,coslats)
indpac_field30_1_ave = global_ave(indpac_field30_1,coslats)

print('temperature')
print(control_temp_ave.mean())
print(pacnino_temp_ave.mean())
print(pacnina_temp_ave.mean())
print(indpiod_temp_ave.mean())
print(indniod_temp_ave.mean())
print(indpac_temp_ave.mean())
print('RH')
print(control_rh_ave.mean())
print(pacnino_rh_ave.mean())
print(pacnina_rh_ave.mean())
print(indpiod_rh_ave.mean())
print(indniod_rh_ave.mean())
print(indpac_rh_ave.mean())
print('cld')
print(control_field30_1_ave.mean())
print(pacnino_field30_1_ave.mean())
print(pacnina_field30_1_ave.mean())
print(indpiod_field30_1_ave.mean())
print(indniod_field30_1_ave.mean())
print(indpac_field30_1_ave.mean())
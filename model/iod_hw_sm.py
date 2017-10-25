# -*- coding: utf-8 -*-
"""
Created on Mon Oct 09 14:40:11 2017

@author: tammas
"""

import numpy as np
import scipy.stats as stats
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
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
# load soil moisture climatology from control
controlnc = nc.MFDataset(control)
lats = controlnc.variables['latitude'][:]
lons = controlnc.variables['longitude'][:]
times = controlnc.variables['t'][:]
dates = nc.num2date(times,units=controlnc.variables['t'].units)
date_range = pd.period_range(dates[0],dates[-1], freq='M')
summer = (date_range.month==12)|(date_range.month==1)|(date_range.month==2)
# Soil moisture
sm_hold = np.squeeze(controlnc.variables['sm'][:,0:3,...].sum(axis=1))
sm_hold = sm_hold[summer]
sm_hold = sm_hold[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
ctrl_sm = np.zeros((30,)+sm_hold.shape[1:])
for yr in xrange(0,30):
    ctrl_sm[yr,...] = sm_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
ctrl_sm = np.ma.array(ctrl_sm,mask=ctrl_sm==0)
sm_summer_clim = np.mean(ctrl_sm, axis=0)


# Load the ensemble soil moisture data
def get_ens(experiment_ens):
    summer_files = glob.glob(experiment_ens[0]+'/*.pa2000-12.nc')
    summernc = nc.MFDataset(summer_files)
    dshape = summernc.variables['sm'].shape[2:]
    summernc.close()
    nens = len(experiment_ens)
    sm = np.zeros((nens,)+dshape)
    for n,ens in enumerate(experiment_ens):
        print ens
        summer_files = glob.glob(ens+'/*.pa2000-12.nc') # December
        summer_files += glob.glob(ens+'/*.pa2001-01.nc') # Janurary
        summer_files += glob.glob(ens+'/*.pa2001-02.nc') # February
        summernc = nc.MFDataset(summer_files)
        sm[n,...] = np.squeeze(np.mean(summernc.variables['sm'][:,0:3,...].sum(axis=1),axis=0))
        summernc.close()
    sm = sm[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
    sm = np.ma.array(sm,mask=sm==0)
    return sm

pacnino_sm = get_ens(pacnino_ens)
pacnina_sm = get_ens(pacnina_ens)
indpiod_sm = get_ens(indpiod_ens)
indniod_sm = get_ens(indniod_ens)
indpac_sm = get_ens(indpac_ens)

lats = lats[(lats<-10.)&(lats>-44.)]
lons = lons[(lons<156.)&(lons>112.)]

# Load heatwave frequency
# Load an ensemble heatwave data.
def load_ensemble_hw(filename, hwdefinition='EHF', get_latlon=False):
    """Load the Australian hw ensemble data and lats and lons.

    Arguments:
    filename -- the path and filename of the file to load.
    hwdefinition -- the heatwave definition the file contains.

    Returns:
    hwf -- frequency
    lats -- latitudes
    lons -- longitudes
    """
    ncfile = nc.Dataset(filename)
    
    lats = ncfile.variables['lat'][:]
    lons = ncfile.variables['lon'][:]

    hwf = ncfile.variables['HWF_'+hwdefinition][0]
    hwf = hwf[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    lats = lats[(lats<-10.)&(lats>-44.)]
    lons = lons[(lons<156.)&(lons>112.)]

    if get_latlon==False:
        lats, lons = 0, 0

    ncfile.close()
    return hwf, lats, lons

# The list of ensebles. vaowf is missing because simulation failed.
control_ensembles = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh', 'vaowf',
        'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp','vaowq','vaowr','vaows',
        'vaowt','vaowu','vaowv','vaoww','vaowx','vaowy','vaowz','vaqgi','vaqgj','vaqgk']
pacnino_ensembles = ['varma','varmc','varme','varmg','varmi','varmk','varmm',
               'varmo','varmq','varms','varmu','varmw','varmy','varna',
               'varnc','varmb','varmd','varmf','varmh','varmj','varml',
               'varmn','varmp','varmr','varmt','varmv','varmx','varmz',
               'varnb','varnd']
pacnina_ensembles = ['varoa','varoc','varoe','varog','varoi','varok','varom',
               'varoo','varoq','varos','varou','varow','varoy','varpa',
               'varpc','varob','varod','varof','varoh','varoj','varol',
               'varon','varop','varor','varot','varov','varox','varoz',
               'varpb','varpd']
indpiod_ensembles = ['varqb','varqd','varqf','varqh','varqj','varql','varqn',
               'varqp','varqr','varqt','varqv','varqx','varqz','varrc',
               'varqa','varqc','varqe','varqg','varqi','varqk','varqm',
               'varqo','varqq','varqs','varqu','varqw','varqy','varrb',
               'varrd']
indniod_ensembles = ['varsa','varsb','varsc','varsd','varse','varsf','varsg',
               'varsh','varsi','varsj','varsk','varsl','varsm','varsn',
               'varso','varsp','varsq','varsr','varss','varst','varsu',
               'varsv','varsx','varsz','varta','vartb','vartc','vartd']
indpac_ensembles = ['vasba','vasbc','vasbe','vasbg','vasbi','vasbk','vasbm',
                  'vasbo','vasbq','vasbs','vasbu','vasca',
                  'vascc','vasbb','vasbd','vasbf','vasbh','vasbj','vasbl',
                  'vasbn','vasbp','vasbr','vasbt','vasbv','vasbx','vasbz',
                  'vascb','vascd']

# Directory of data
directory = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
filename = directory+'vamrk/EHF_heatwaves_ACCESS1.3_vamrk_yearly_summer.nc'

# Get lats and lons
_, lats, lons = load_ensemble_hw(filename, get_latlon=True)

# Initialise arrays
empty = np.ma.ones((len(pacnino_ensembles), len(lats), len(lons)))*np.nan
pacnino_hwf = empty.copy()
empty = np.ma.ones((len(pacnina_ensembles), len(lats), len(lons)))*np.nan
pacnina_hwf = empty.copy()
empty = np.ma.ones((len(indpiod_ensembles), len(lats), len(lons)))*np.nan
indpiod_hwf = empty.copy()
empty = np.ma.ones((len(indniod_ensembles), len(lats), len(lons)))*np.nan
indniod_hwf = empty.copy()
empty = np.ma.ones((len(indpac_ensembles), len(lats), len(lons)))*np.nan
indpac_hwf = empty.copy()

# Load control heatwave data
filename = directory+'vamrb'+'/EHF_heatwaves_ACCESS1.3_vamrb_yearly_summer.nc'
controlhwnc = nc.Dataset(filename)
ctrl_hwf = controlhwnc.variables['HWF_EHF'][-30:,...]
lats = controlhwnc.variables['lat'][:]
lons = controlhwnc.variables['lon'][:]
ctrl_hwf = ctrl_hwf[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
controlhwnc.close()
lats = lats[(lats<-10.)&(lats>-44.)]
lons = lons[(lons<156.)&(lons>112.)]
directory = '/srv/ccrc/data48/z5032520/ehfheatwaves/'
# Load pacnino ensebles
for n, ens in enumerate(pacnino_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    pacnino_hwf[n], _, _ = load_ensemble_hw(filename)
# Load pacnina ensebles
for n, ens in enumerate(pacnina_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    pacnina_hwf[n], _, _ = load_ensemble_hw(filename)
# Load indpiod ensebles
for n, ens in enumerate(indpiod_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indpiod_hwf[n], _, _ = load_ensemble_hw(filename)
# Load indniod ensebles
for n, ens in enumerate(indniod_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indniod_hwf[n], _, _ = load_ensemble_hw(filename)
# Load indpac ensebles
for n, ens in enumerate(indpac_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indpac_hwf[n], _, _ = load_ensemble_hw(filename)

print "Done Loading"

# Define the regions of interest 
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)
swaus = (115,-27)
seaus = (141.,-31.)
neaus = (139.,-18.)
naus = (129.,-12)
eaus = (145.5,-24)

# Cut out the regions of interest
def cut_region(data,region,lats,lons):
    size = 7.5
    xi = (lons>region[0])&(lons<region[0]+size)
    yi = (lats<region[1])&(lats>region[1]-size)
    return data[:,yi,:][:,:,xi].mean(axis=2).mean(axis=1)

ctrl_neaus_hwf = cut_region(ctrl_hwf,neaus,lats,lons)
ctrl_seaus_hwf = cut_region(ctrl_hwf,seaus,lats,lons) 
ctrl_naus_hwf = cut_region(ctrl_hwf,naus,lats,lons)
ctrl_eaus_hwf = cut_region(ctrl_hwf,eaus,lats,lons)
ctrl_swaus_hwf = cut_region(ctrl_hwf,swaus,lats,lons)

pacnino_neaus_hwf = cut_region(pacnino_hwf,neaus,lats,lons)
pacnino_seaus_hwf = cut_region(pacnino_hwf,seaus,lats,lons) 
pacnino_naus_hwf = cut_region(pacnino_hwf,naus,lats,lons)
pacnino_eaus_hwf = cut_region(pacnino_hwf,eaus,lats,lons)
pacnino_swaus_hwf = cut_region(pacnino_hwf,swaus,lats,lons)

pacnina_neaus_hwf = cut_region(pacnina_hwf,neaus,lats,lons)
pacnina_seaus_hwf = cut_region(pacnina_hwf,seaus,lats,lons) 
pacnina_naus_hwf = cut_region(pacnina_hwf,naus,lats,lons)
pacnina_eaus_hwf = cut_region(pacnina_hwf,eaus,lats,lons)
pacnina_swaus_hwf = cut_region(pacnina_hwf,swaus,lats,lons)

indpiod_neaus_hwf = cut_region(indpiod_hwf,neaus,lats,lons)
indpiod_seaus_hwf = cut_region(indpiod_hwf,seaus,lats,lons) 
indpiod_naus_hwf = cut_region(indpiod_hwf,naus,lats,lons)
indpiod_eaus_hwf = cut_region(indpiod_hwf,eaus,lats,lons)
indpiod_swaus_hwf = cut_region(indpiod_hwf,swaus,lats,lons)

indniod_neaus_hwf = cut_region(indniod_hwf,neaus,lats,lons)
indniod_seaus_hwf = cut_region(indniod_hwf,seaus,lats,lons) 
indniod_naus_hwf = cut_region(indniod_hwf,naus,lats,lons)
indniod_eaus_hwf = cut_region(indniod_hwf,eaus,lats,lons)
indniod_swaus_hwf = cut_region(indniod_hwf,swaus,lats,lons)

indpac_neaus_hwf = cut_region(indpac_hwf,neaus,lats,lons)
indpac_seaus_hwf = cut_region(indpac_hwf,seaus,lats,lons) 
indpac_naus_hwf = cut_region(indpac_hwf,naus,lats,lons)
indpac_eaus_hwf = cut_region(indpac_hwf,eaus,lats,lons)
indpac_swaus_hwf = cut_region(indpac_hwf,swaus,lats,lons)

ctrl_neaus_sm = cut_region(ctrl_sm,neaus,lats,lons)
ctrl_seaus_sm = cut_region(ctrl_sm,seaus,lats,lons) 
ctrl_naus_sm = cut_region(ctrl_sm,naus,lats,lons)
ctrl_eaus_sm = cut_region(ctrl_sm,eaus,lats,lons)
ctrl_swaus_sm = cut_region(ctrl_sm,swaus,lats,lons)

pacnino_neaus_sm = cut_region(pacnino_sm,neaus,lats,lons)
pacnino_seaus_sm = cut_region(pacnino_sm,seaus,lats,lons) 
pacnino_naus_sm = cut_region(pacnino_sm,naus,lats,lons)
pacnino_eaus_sm = cut_region(pacnino_sm,eaus,lats,lons)
pacnino_swaus_sm = cut_region(pacnino_sm,swaus,lats,lons)

pacnina_neaus_sm = cut_region(pacnina_sm,neaus,lats,lons)
pacnina_seaus_sm = cut_region(pacnina_sm,seaus,lats,lons) 
pacnina_naus_sm = cut_region(pacnina_sm,naus,lats,lons)
pacnina_eaus_sm = cut_region(pacnina_sm,eaus,lats,lons)
pacnina_swaus_sm = cut_region(pacnina_sm,swaus,lats,lons)

indpiod_neaus_sm = cut_region(indpiod_sm,neaus,lats,lons)
indpiod_seaus_sm = cut_region(indpiod_sm,seaus,lats,lons) 
indpiod_naus_sm = cut_region(indpiod_sm,naus,lats,lons)
indpiod_eaus_sm = cut_region(indpiod_sm,eaus,lats,lons)
indpiod_swaus_sm = cut_region(indpiod_sm,swaus,lats,lons)

indniod_neaus_sm = cut_region(indniod_sm,neaus,lats,lons)
indniod_seaus_sm = cut_region(indniod_sm,seaus,lats,lons) 
indniod_naus_sm = cut_region(indniod_sm,naus,lats,lons)
indniod_eaus_sm = cut_region(indniod_sm,eaus,lats,lons)
indniod_swaus_sm = cut_region(indniod_sm,swaus,lats,lons)

indpac_neaus_sm = cut_region(indpac_sm,neaus,lats,lons)
indpac_seaus_sm = cut_region(indpac_sm,seaus,lats,lons) 
indpac_naus_sm = cut_region(indpac_sm,naus,lats,lons)
indpac_eaus_sm = cut_region(indpac_sm,eaus,lats,lons)
indpac_swaus_sm = cut_region(indpac_sm,swaus,lats,lons)

# Plot the data
#sl,inter,rv,pv,std = stats.theilslopes(ctrl_neaus_sm,ctrl_neaus_hwf)

plt.figure()
plt.scatter(ctrl_neaus_sm,ctrl_neaus_hwf,c='k',marker='.',label='Control')
plt.scatter(pacnino_neaus_sm,pacnino_neaus_hwf,c='r',marker='+',label='El Nino')
plt.scatter(pacnina_neaus_sm,pacnina_neaus_hwf,c='b',marker='+',label='La Nina')
plt.scatter(indpiod_neaus_sm,indpiod_neaus_hwf,c='m',marker='^',label='+IOD')
plt.scatter(indniod_neaus_sm,indniod_neaus_hwf,c='c',marker='^',label='-IOD')
plt.scatter(indpac_neaus_sm,indpac_neaus_hwf,c='brown',marker='o',label='El Nino +IOD')
#sl,inter,_,_ = stats.theilslopes(ctrl_neaus_hwf,ctrl_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'k')
#sl,inter,_,_ = stats.theilslopes(pacnino_neaus_hwf,pacnino_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'r')
#sl,inter,_,_ = stats.theilslopes(pacnina_neaus_hwf,pacnina_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'b')
#sl,inter,_,_ = stats.theilslopes(indpiod_neaus_hwf,indpiod_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'m')
#sl,inter,_,_ = stats.theilslopes(indniod_neaus_hwf,indniod_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'c')
#sl,inter,_,_ = stats.theilslopes(indpac_neaus_hwf,indpac_neaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'brown')
plt.xlabel('Soil moisture ($kg/m^{2}$)')
plt.ylabel('HWF (days)')
plt.title('Northeast Australia')
plt.legend()
plt.savefig('neaus_hwf_sm.png',format='png',dpi=250)

plt.figure()
plt.scatter(ctrl_naus_sm,ctrl_naus_hwf,c='k',marker='.',label='Control')
plt.scatter(pacnino_naus_sm,pacnino_naus_hwf,c='r',marker='+',label='El Nino')
plt.scatter(pacnina_naus_sm,pacnina_naus_hwf,c='b',marker='+',label='La Nina')
plt.scatter(indpiod_naus_sm,indpiod_naus_hwf,c='m',marker='^',label='+IOD')
plt.scatter(indniod_naus_sm,indniod_naus_hwf,c='c',marker='^',label='-IOD')
plt.scatter(indpac_naus_sm,indpac_naus_hwf,c='brown',marker='o',label='El Nino +IOD')
#sl,inter,_,_ = stats.theilslopes(ctrl_naus_hwf,ctrl_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'k')
#sl,inter,_,_ = stats.theilslopes(pacnino_naus_hwf,pacnino_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'r')
#sl,inter,_,_ = stats.theilslopes(pacnina_naus_hwf,pacnina_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'b')
#sl,inter,_,_ = stats.theilslopes(indpiod_naus_hwf,indpiod_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'m')
#sl,inter,_,_ = stats.theilslopes(indniod_naus_hwf,indniod_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'c')
#sl,inter,_,_ = stats.theilslopes(indpac_naus_hwf,indpac_naus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'brown')
plt.xlabel('Soil moisture ($kg/m^{2}$)')
plt.ylabel('HWF (days)')
plt.title('North Australia')
plt.legend()
plt.savefig('naus_hwf_sm.png',format='png',dpi=250)

plt.figure()
plt.scatter(ctrl_eaus_sm,ctrl_eaus_hwf,c='k',marker='.',label='Control')
plt.scatter(pacnino_eaus_sm,pacnino_eaus_hwf,c='r',marker='+',label='El Nino')
plt.scatter(pacnina_eaus_sm,pacnina_eaus_hwf,c='b',marker='+',label='La Nina')
plt.scatter(indpiod_eaus_sm,indpiod_eaus_hwf,c='m',marker='^',label='+IOD')
plt.scatter(indniod_eaus_sm,indniod_eaus_hwf,c='c',marker='^',label='-IOD')
plt.scatter(indpac_eaus_sm,indpac_eaus_hwf,c='brown',marker='o',label='El Nino +IOD')
#sl,inter,_,_ = stats.theilslopes(ctrl_eaus_hwf,ctrl_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'k')
#sl,inter,_,_ = stats.theilslopes(pacnino_eaus_hwf,pacnino_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'r')
#sl,inter,_,_ = stats.theilslopes(pacnina_eaus_hwf,pacnina_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'b')
#sl,inter,_,_ = stats.theilslopes(indpiod_eaus_hwf,indpiod_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'m')
#sl,inter,_,_ = stats.theilslopes(indniod_eaus_hwf,indniod_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'c')
#sl,inter,_,_ = stats.theilslopes(indpac_eaus_hwf,indpac_eaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'brown')
plt.xlabel('Soil moisture ($kg/m^{2}$)')
plt.ylabel('HWF (days)')
plt.title('East Australia')
plt.legend()
plt.savefig('eaus_hwf_sm.png',format='png',dpi=250)

plt.figure()
plt.scatter(ctrl_seaus_sm,ctrl_seaus_hwf,c='k',marker='.',label='Control')
plt.scatter(pacnino_seaus_sm,pacnino_seaus_hwf,c='r',marker='+',label='El Nino')
plt.scatter(pacnina_seaus_sm,pacnina_seaus_hwf,c='b',marker='+',label='La Nina')
plt.scatter(indpiod_seaus_sm,indpiod_seaus_hwf,c='m',marker='^',label='+IOD')
plt.scatter(indniod_seaus_sm,indniod_seaus_hwf,c='c',marker='^',label='-IOD')
plt.scatter(indpac_seaus_sm,indpac_seaus_hwf,c='brown',marker='o',label='El Nino +IOD')
#sl,inter,_,_ = stats.theilslopes(ctrl_seaus_hwf,ctrl_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'k')
#sl,inter,_,_ = stats.theilslopes(pacnino_seaus_hwf,pacnino_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'r')
#sl,inter,_,_ = stats.theilslopes(pacnina_seaus_hwf,pacnina_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'b')
#sl,inter,_,_ = stats.theilslopes(indpiod_seaus_hwf,indpiod_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'m')
#sl,inter,_,_ = stats.theilslopes(indniod_seaus_hwf,indniod_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'c')
#sl,inter,_,_ = stats.theilslopes(indpac_seaus_hwf,indpac_seaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'brown')
plt.xlabel('Soil moisture ($kg/m^{2}$)')
plt.ylabel('HWF (days)')
plt.title('Southeast Australia')
plt.legend()
plt.savefig('seaus_hwf_sm.png',format='png',dpi=250)

plt.figure()
plt.scatter(ctrl_swaus_sm,ctrl_seaus_hwf,c='k',marker='.',label='Control')
plt.scatter(pacnino_swaus_sm,pacnino_swaus_hwf,c='r',marker='+',label='El Nino')
plt.scatter(pacnina_swaus_sm,pacnina_swaus_hwf,c='b',marker='+',label='La Nina')
plt.scatter(indpiod_swaus_sm,indpiod_swaus_hwf,c='m',marker='^',label='+IOD')
plt.scatter(indniod_swaus_sm,indniod_swaus_hwf,c='c',marker='^',label='-IOD')
plt.scatter(indpac_swaus_sm,indpac_swaus_hwf,c='brown',marker='o',label='El Nino +IOD')
#sl,inter,_,_ = stats.theilslopes(ctrl_swaus_hwf,ctrl_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'k')
#sl,inter,_,_ = stats.theilslopes(pacnino_swaus_hwf,pacnino_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'r')
#sl,inter,_,_ = stats.theilslopes(pacnina_swaus_hwf,pacnina_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'b')
#sl,inter,_,_ = stats.theilslopes(indpiod_swaus_hwf,indpiod_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'m')
#sl,inter,_,_ = stats.theilslopes(indniod_swaus_hwf,indniod_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'c')
#sl,inter,_,_ = stats.theilslopes(indpac_swaus_hwf,indpac_swaus_sm)
#plt.plot([1.5,4.5],[sl*1.5+inter,sl*4.5+inter],'brown')
plt.xlabel('Soil moisture ($kg/m^{2}$)')
plt.ylabel('HWF (days)')
plt.title('Southwest Australia')
plt.legend()
plt.savefig('swaus_hwf_sm.png',format='png',dpi=250)
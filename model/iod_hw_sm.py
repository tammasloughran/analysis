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

# load soil moisture climatology from control
controlnc = nc.MFDataset(control)
lats = controlnc.variables['latitude'][:]
lons = controlnc.variables['longitude'][:]
times = controlnc.variables['t'][:]
dates = nc.num2date(times,units=controlnc.variables['t'].units)
date_range = pd.period_range(dates[0],dates[-1], freq='M')
summer = (date_range.month==12)|(date_range.month==1)|(date_range.month==2)
# Soil moisture
sm_hold = np.squeeze(controlnc.variables['sm'][:,5,...])
sm_hold = sm_hold[summer]
sm_hold = sm_hold[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
sm_ctrl = np.zeros((30,)+sm_hold.shape[1:])
for yr in xrange(0,30):
    sm_ctrl[yr,...] = sm_hold[(yr*3)+2:(yr*3)+5,...].mean(axis=0)
sm_summer_clim = np.mean(sm_ctrl, axis=0)


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
        sm[n,...] = np.squeeze(np.mean(summernc.variables['sm'][:,5,...],axis=0))
        summernc.close()
    sm = sm[:,(lats<-10.)&(lats>-44.),:][:,:,(lons<156.)&(lons>112.)]
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
control_hwf = controlhwnc.variables['HWF_EHF'][-30:,...]
controlhwnc.close()
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



# Define the regions of interest 
#(129.,-12),(139.,-18),(145.5,-24),(141,-31),(115,-27)
seaus = (141.,-31.)
neaus = (139.,-18.)
naus = (129.,-12)
eaus = (145.5,-24)



print "Done Loading"
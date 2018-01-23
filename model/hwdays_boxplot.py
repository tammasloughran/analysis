# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:42:05 2018

@author: tammas
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


# Load an ensemble heatwave data.
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


# The list of ensebles. vaowf is missing because simulation failed.
control_ensembles = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh', 'vaowf',
        'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp','vaowq','vaowr','vaows',
        'vaowt','vaowu','vaowv','vaoww','vaowx','vaowy','vaowz','vaqgi','vaqgj','vaqgk']
nino_ensembles = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg', 'vaoqh',
        'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq','vaoqr', 'vaoqn',
        'vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy','vaoqz','vaqgl','vaqgm','vaqgn']
nina_ensembles = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl','vamrm',
        'vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt','vamru','vamrv','vamrw','vamrx',
        'vamry','vamrz','vaqga','vaqgb','vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']
modoki_ensembles = ['vaqoc','vaqog','vaqok','vaqoo','vaqos','vaqow','vaqpa','vaqod',
              'vaqoh','vaqol','vaqop','vaqot','vaqox','vaqpb','vaqoa','vaqoe',
              'vaqoi','vaqom','vaqoq','vaqou','vaqoy','vaqpc','vaqob','vaqof',
              'vaqoj','vaqon','vaqor','vaqov','vaqoz','vaqpd']
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
_,_,_,_,_,_, lats, lons = load_ensemble_hw(filename, get_latlon=True)

# Initialise arrays
empty = np.ma.ones((len(control_ensembles), len(lats), len(lons)))*np.nan
control = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(nino_ensembles), len(lats), len(lons)))*np.nan
elnino = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(nina_ensembles), len(lats), len(lons)))*np.nan
lanina = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(modoki_ensembles), len(lats), len(lons)))*np.nan
modoki = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(pacnino_ensembles), len(lats), len(lons)))*np.nan
pacnino = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(pacnina_ensembles), len(lats), len(lons)))*np.nan
pacnina = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(indpiod_ensembles), len(lats), len(lons)))*np.nan
indpiod = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(indniod_ensembles), len(lats), len(lons)))*np.nan
indniod = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}
empty = np.ma.ones((len(indpac_ensembles), len(lats), len(lons)))*np.nan
indpac = {'hwf': empty.copy(), 'hwn': empty.copy(), \
        'hwd': empty.copy(), 'hwa': empty.copy(), \
        'hwm': empty.copy(), 'hwt': empty.copy()}

# Load control ensebles
for n, ens in enumerate(control_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    control['hwf'][n], control['hwn'][n], \
    control['hwd'][n], control['hwa'][n], \
    control['hwm'][n], control['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load el nino ensebles
for n, ens in enumerate(nino_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    elnino['hwf'][n], elnino['hwn'][n], \
    elnino['hwd'][n], elnino['hwa'][n], \
    elnino['hwm'][n], elnino['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load la nina ensebles
for n, ens in enumerate(nina_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    lanina['hwf'][n], lanina['hwn'][n], \
    lanina['hwd'][n], lanina['hwa'][n], \
    lanina['hwm'][n], lanina['hwt'][n], _, _ = load_ensemble_hw(filename)

# Regions
seaus = (141.,-31.)
neaus = (139.,-18.)
eaus = (145.5,-24)
naus = (129.,-12)

def index_region(x,y,data):
    ys = (lats>=y-7.5)&(lats<=y)
    xs = (lons<=x+7.5)&(lons>=x)
    return data[:,ys,:][:,:,xs].mean(axis=1).mean(axis=1)

# Control neaus climatology
neaushwf_ctrl = index_region(neaus[0],neaus[1],control['hwf'])

seaushwf_nino = index_region(seaus[0],seaus[1],elnino['hwf'])
eaushwf_nino = index_region(eaus[0],eaus[1],elnino['hwf'])
neaushwf_nino = index_region(neaus[0],neaus[1],elnino['hwf'])
naushwf_nino = index_region(naus[0],naus[1],elnino['hwf'])

seaushwf_nina = index_region(seaus[0],seaus[1],lanina['hwf'])
eaushwf_nina = index_region(eaus[0],eaus[1],lanina['hwf'])
neaushwf_nina = index_region(neaus[0],neaus[1],lanina['hwf'])
naushwf_nina = index_region(naus[0],naus[1],lanina['hwf'])

ninodata = [[],naushwf_nino,neaushwf_nino,eaushwf_nino,seaushwf_nino]
ninadata = [naushwf_nina,neaushwf_nina,eaushwf_nina,seaushwf_nina]

ninopos = [1,2,5,8,11] 
ninapos = [1,4,7,10]
lbls = ['North','Northeast','East','Southeast']

nino_box = plt.boxplot(ninodata,positions=ninopos,whis=100)
nina_box = plt.boxplot(ninadata,positions=ninapos,labels=lbls,whis=100)

for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(nina_box[element], color='b')
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(nino_box[element], color='r')

for x,data in zip(ninopos,ninodata):
    if data==[]: continue
    plt.scatter(np.ones(len(data))*x, data, color='r', marker='+', linewidth=0.5)
for x,data in zip(ninapos,ninadata):
    if data==[]: continue
    plt.scatter(np.ones(len(data))*x, data, color='b', marker='+', linewidth=0.5)

plt.scatter(2,neaushwf_ctrl.mean(), color='k', marker='o')

plt.ylabel('No. heatwave days')
plt.xticks([1.5,4.5,7.5,10.5])
plt.xlim(0,12)
plt.show()




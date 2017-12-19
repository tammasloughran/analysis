# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:52:06 2017

@author: tammas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import scipy.stats as stats

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
directory = '/srv/ccrc/data48/z5032520/ehfheatwaves/'
# Load pacnino ensebles
for n, ens in enumerate(pacnino_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    pacnino['hwf'][n], pacnino['hwn'][n], \
    pacnino['hwd'][n], pacnino['hwa'][n], \
    pacnino['hwm'][n], pacnino['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load pacnina ensebles
for n, ens in enumerate(pacnina_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    pacnina['hwf'][n], pacnina['hwn'][n], \
    pacnina['hwd'][n], pacnina['hwa'][n], \
    pacnina['hwm'][n], pacnina['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load indpiod ensebles
for n, ens in enumerate(indpiod_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indpiod['hwf'][n], indpiod['hwn'][n], \
    indpiod['hwd'][n], indpiod['hwa'][n], \
    indpiod['hwm'][n], indpiod['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load indniod ensebles
for n, ens in enumerate(indniod_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indniod['hwf'][n], indniod['hwn'][n], \
    indniod['hwd'][n], indniod['hwa'][n], \
    indniod['hwm'][n], indniod['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load indpac ensebles
for n, ens in enumerate(indpac_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    indpac['hwf'][n], indpac['hwn'][n], \
    indpac['hwd'][n], indpac['hwa'][n], \
    indpac['hwm'][n], indpac['hwt'][n], _, _ = load_ensemble_hw(filename)

# Load HWF EOF1
#nea = np.load('mkdnino_HWF_masks.npy')
#nea[0,24:26,:] = 0
#nea = nea[0].astype('bool')


def plotmap_test(data,aspect,test,name,colours='RdBu_r'):
    #fmt = '%1.0f'
    plt.figure()
    parallels = np.arange(-40., -9., 10.)
    meridians = np.arange(120., 160., 10.,)
    map_axes = Basemap(projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    xx, yy = np.meshgrid(lons, lats)
    x, y = map_axes(xx,yy)
    xx = xx - 0.5 # The data projection is slightly off compared to the coastlines.
    yy = yy - 0.5
    px, py = map_axes(xx,yy)
    data = np.ma.array(data,mask=np.isnan(data))
    map_axes.pcolormesh(px,py,data,cmap=colours,vmin=0,vmax=1)
    map_axes.drawparallels(parallels, labels=[True,False,False,False], linewidth=0.5)
    map_axes.drawmeridians(meridians, labels=[False,False,False,True], linewidth=0.5)
    map_axes.drawcoastlines()
    plt.title(test+' '+aspect)
    plt.savefig(name, format='eps')
    plt.close()

# Load the mask
msk = control['hwf'].mask[0,...]

def plotmaps(data,signif,name,filename,colours='viridis'):
    fmt = '%1.0f'
    if name=='hwf':
        units = 'Days'
        cints = np.arange(-15,16,3)
    elif name=='hwd':
        units = 'Days'
        cints = np.arange(-4,5,1)
    elif name=='hwn':
        units = 'Heatwaves'
        cints = np.arange(-2,2.4,0.4)
        fmt = '%1.1f'
    elif name=='hwa':
        units = '$^{\circ}C^{2}$'
        cints = np.arange(-15,16,3)
    elif name=='hwm':
        units = '$^{\circ}C^{2}$'
        cints = np.arange(-5,6,1)
    elif name=='hwt':
        units = 'Days'
        cints = np.arange(-30,31,10)
    else: 
        units = 'stuff'
    plt.figure()
    parallels = np.arange(-40., -9., 10.)
    meridians = np.arange(120., 160., 10.,)
    map_axes = Basemap(projection='cyl',
        llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156,resolution='l')
    map_axes.drawparallels(parallels, labels=[True,False,False,False])
    map_axes.drawmeridians(meridians, labels=[False,False,False,True])
    xx, yy = np.meshgrid(lons, lats)
    x, y = map_axes(xx,yy)
    xx = xx - 0.5 # The data projection is slightly off compared to the coastlines.
    yy = yy - 0.5
    px, py = map_axes(xx,yy)
    data = np.ma.array(data,mask=np.isnan(data))
    shade = map_axes.pcolormesh(px,py,data,cmap=colours,vmin=cints[0],vmax=cints[-1])
    map_axes.contourf(x,y,signif, 1, colors='none', hatches=[None,'x'])
    cont = map_axes.contour(x,y,data,levels=cints,colors='k',linewidths=0.5)
    for c in cont.collections:
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    #plt.clabel(cont, fmt=fmt, fontsize=10)
    cb = map_axes.colorbar(shade, location='bottom', pad=0.25, ticks=cints)
    cb.ax.set_xlabel(units)
    map_axes.drawparallels(parallels, labels=[True,False,False,False])
    map_axes.drawmeridians(meridians, labels=[False,False,False,True])
    map_axes.drawcoastlines()
    plt.title(name)
    plt.tight_layout()
    plt.savefig(filename, format='png', dpi=200)
    plt.close()


def kruskal_3d(ad,bd):
    """Does a Kruskal test for 3 dimensional data allong the first axis.
    ad and bd are the two groups to test. 
    Returns hs, the H statistic and pv the p-value for a ? sided hypothesis test.
    """
    hs = np.ones(msk.shape)*np.nan
    hs = np.ma.array(hs,mask=msk)
    pv = hs.copy()
    for y in xrange(ad.shape[1]):
        for x in xrange(bd.shape[2]):
            if hs.mask[y,x]==False:
                hs[y,x], pv[y,x] = stats.kruskal(ad[:,y,x],bd[:,y,x],nan_policy='omit')
    return hs, pv

def ks_3d(ad,bd):
    """Does a KS test for 3 dimensional data allong the first axis.
    ad and bd are the two groups to test. 
    Returns ks, the KS statistic and pv the p-value for a 2 sided hypothesis test.
    """
    ks = np.ones(msk.shape)*np.nan
    ks = np.ma.array(ks,mask=msk)
    pv = ks.copy()
    for y in xrange(ad.shape[1]):
        for x in xrange(bd.shape[2]):
            if ks.mask[y,x]==False:
                ks[y,x], pv[y,x] = stats.ks_2samp(ad[:,y,x],bd[:,y,x])
    return ks, pv

def mannwhitneyu_3d(ad,bd):
    """Does a Mann-Whitney U test for 3 dimensional data allong the first axis.
    ad and bd are the two groups to test. 
    Returns us, the U statistic and pv the p-value for a two sided hypothesis test.
    """
    us = np.ones(msk.shape)*np.nan
    us = np.ma.array(us,mask=msk)
    pv = us.copy()
    for y in xrange(ad.shape[1]):
        for x in xrange(bd.shape[2]):
            if us.mask[y,x]==False:
                us[y,x], pv[y,x] = stats.mannwhitneyu(ad[:,y,x],bd[:,y,x], alternative='two-sided')
    return us, pv


#def random_replace_subsample(data, nsub=20, average=True):
#    """Select a random subsample of size nsub from data and optionally 
#    calculate the average.
#    
#    Inputs
#    data - Data to subsample
#    nsub - the number of elements in the subsample
#    average - Returns the average allong the first axis if true
#    
#    Returns
#    sample - the subsampled data
#    """
#    sample = np.ma.zeros((nsub,)+data.shape[1:])
#    for n in xrange(nsub):
#        sample[n] = data[int(np.random.rand()*data.shape[0])]
#    if average==True: sample = sample.mean(axis=0)
#    return sample
#
#
#def make_pdf(data, obs, title=''):
#    distribution = data[:,nea].mean(axis=1)
#    plt.figure()
#    plt.hist(distribution, bins=50, normed=True, color='grey')
#    plt.xticks(np.arange(-20,20,2))
#    plt.vlines(distribution.mean(), 0, 0.25, colors='k')
#    plt.vlines(np.percentile(distribution, 5), 0, 0.25, colors='k', linestyles='dashed')
#    plt.vlines(np.percentile(distribution, 95), 0, 0.25, colors='k', linestyles='dashed')
#    plt.vlines(obs, 0, 0.25, colors='k', linestyles='dotted')
#    plt.xlabel('Days')
#    ylbl = plt.ylabel('Pr.')
#    ylbl.set_rotation(0)
#    plt.title(title)
#    plt.show()
#
#
#def bootstrap(data1, data2, obs, nsamples=10000, title=''):
#    """Calculate the fifth and ninetyfifth percentil confidence intervals 
#    for the difference of means of two groups using bootstrapping with 
#    replacement. 
#    
#    Inputs
#    data1, data2 - the two groups to calculate the difference of means
#    nsample - the number of bootstrap samples to perform
#    
#    Returns
#    avg - the average difference in the means
#    fifthint - the fifth percentile confidence interval
#    ninfifthint - the ninety-fifth percentile confidence interval
#    """
#    dffrnc = np.ma.zeros((nsamples,)+data1.shape[1:])
#    for ni in xrange(nsamples):
#        dffrnc[ni,...] = random_replace_subsample(data1, nsub=15) - random_replace_subsample(data2, nsub=15)
#    avg = dffrnc.mean(axis=0)
#    make_pdf(dffrnc,obs,title=title)
#    fifthint = np.percentile(dffrnc, 5, axis=0)
#    ninfifthint = np.percentile(dffrnc, 95, axis=0)
#    return avg, fifthint, ninfifthint
    
    
## Figure for chapter 2
#def subplot_aus(data, lons, lats, signif, ax):
#    map_axes = Basemap(ax=ax, projection='cyl',
#        llcrnrlat=-44, urcrnrlat=-10,
#        llcrnrlon=112, urcrnrlon=156,resolution='l')
#    xx, yy = np.meshgrid(lons, lats)
#    x, y = map_axes(xx,yy)
#    xx = xx - 0.5 # The data projection is slightly off compared to the coastlines.
#    yy = yy - 0.5
#    px, py = map_axes(xx,yy)
#    data = np.ma.array(data,mask=np.isnan(data))
#    shade = map_axes.pcolormesh(px,py,data,cmap='RdBu_r',vmin=cints[0],vmax=cints[-1])
#    map_axes.contourf(x,y,signif, 1, colors='none', hatches=[None,'xx'])
#    cont = map_axes.contour(x,y,data,colors='k',linewidths=0.5,levels=cints)
#    for c in cont.collections:
#        if c.get_linestyle() == [(None, None)]:
#           continue
#        else:
#            c.set_dashes([(0, (2.0, 2.0))])
#    map_axes.drawcoastlines()
#    return map_axes, shade
#    #map_axes.colorbar(shade,location='bottom',ax=axes[2][i])

threshold = 0.05
# Nino diff
for aspect in control:
    print aspect
    # Differene of means for elnino
    # pacnino
    ts, pv = stats.ttest_ind(pacnino[aspect],control[aspect],axis=0,equal_var=False, nan_policy='omit')
    sig = pv<threshold
    pacnino_diff = np.ma.mean(pacnino[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(pacnino_diff, sig, aspect, 'mean_'+aspect+'_pacnino_diff.png', 'RdBu_r')
    # pacnina
    ts, pv = stats.ttest_ind(pacnina[aspect],control[aspect],axis=0,equal_var=False, nan_policy='omit')
    sig = pv<threshold
    pacnina_diff = np.ma.mean(pacnina[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(pacnina_diff, sig, aspect, 'mean_'+aspect+'_pacnina_diff.png', 'RdBu_r')
    # indpiod
    ts, pv = stats.ttest_ind(indpiod[aspect],control[aspect],axis=0,equal_var=False, nan_policy='omit')
    sig = pv<threshold
    indpiod_diff = np.ma.mean(indpiod[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(indpiod_diff, sig, aspect, 'mean_'+aspect+'_indpiod_diff.png', 'RdBu_r')
    # indniod
    ts, pv = stats.ttest_ind(indniod[aspect],control[aspect],axis=0,equal_var=False, nan_policy='omit')
    sig = pv<threshold
    indniod_diff = np.ma.mean(indniod[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(indniod_diff, sig, aspect, 'mean_'+aspect+'_indniod_diff.png', 'RdBu_r')
    # indpac
    ts, pv = stats.ttest_ind(indpac[aspect],control[aspect],axis=0,equal_var=False, nan_policy='omit')
    sig = pv<threshold
    indpac_diff = np.ma.mean(indpac[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(indpac_diff, sig, aspect, 'mean_'+aspect+'_indpac_diff.png', 'RdBu_r')
    
    sig[...] = 0
    plotmaps(pacnino_diff+indpiod_diff,sig,aspect,'ninoppiod_'+aspect+'.png','RdBu_r')


#    # Difference of medians nino
#    # pacnino
#    ad = pacnino[aspect].data
#    bd = control[aspect].data
#    hs, pv = mannwhitneyu_3d(ad,bd)
#    sig = pv<threshold
#    nina_diff = np.ma.median(pacnino[aspect], axis=0)-np.ma.median(control[aspect], axis=0)
#    plotmaps(pacnino_diff, sig, aspect, 'median_'+aspect+'_pacnino_diff.eps', 'RdBu_r')
#    # pacnina
#    ad = pacnina[aspect].data
#    bd = control[aspect].data
#    hs, pv = mannwhitneyu_3d(ad,bd)
#    sig = pv<threshold
#    pacnina_diff = np.ma.median(pacnina[aspect], axis=0)-np.ma.median(control[aspect], axis=0)
#    plotmaps(pacnina_diff, sig, aspect, 'median_'+aspect+'_pacnina_diff.eps', 'RdBu_r')
#    # indpiod
#    ad = indpiod[aspect].data
#    bd = control[aspect].data
#    hs, pv = mannwhitneyu_3d(ad,bd)
#    sig = pv<threshold
#    indpiod_diff = np.ma.median(indpiod[aspect], axis=0)-np.ma.median(control[aspect], axis=0)
#    plotmaps(indpiod_diff, sig, aspect, 'median_'+aspect+'_indpiod_diff.eps', 'RdBu_r')
#    # indniod
#    ad = indniod[aspect].data
#    bd = control[aspect].data
#    hs, pv = mannwhitneyu_3d(ad,bd)
#    sig = pv<threshold
#    indniod_diff = np.ma.median(indniod[aspect], axis=0)-np.ma.median(control[aspect], axis=0)
#    plotmaps(indniod_diff, sig, aspect, 'median_'+aspect+'_indniod_diff.eps', 'RdBu_r')
#    # indpac
#    ad = indpac[aspect].data
#    bd = control[aspect].data
#    hs, pv = mannwhitneyu_3d(ad,bd)
#    sig = pv<threshold
#    indpac_diff = np.ma.median(indpac[aspect], axis=0)-np.ma.median(control[aspect], axis=0)
#    plotmaps(indpac_diff, sig, aspect, 'median_'+aspect+'_indpac_diff.eps', 'RdBu_r')

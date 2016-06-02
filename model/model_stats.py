import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import scipy.stats as stats
import pdb

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

    hwa = ncfile.variables['HWA_'+hwdefinition][0]
    hwa = hwa[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    hwm = ncfile.variables['HWM_'+hwdefinition][0]
    hwm = hwm[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    hwt = ncfile.variables['HWT_'+hwdefinition][0]
    hwt = hwt[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]

    lats = lats[(lats<-10.)&(lats>-44.)]
    lons = lons[(lons<156.)&(lons>112.)]

    if get_latlon==False:
        lats, lons = 0, 0

    return hwf, hwn, hwd, hwa, hwm, hwt, lats, lons

# The list of ensebles. vaowf is missing because simulation failed.
control_ensembles = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh', #'vaowf',
        'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp','vaowq','vaowr','vaows',
        'vaowt','vaowu','vaowv','vaoww','vaowx','vaowy','vaowz','vaqgi','vaqgj','vaqgk']
nino_ensembles = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg', #'vaoqh',
        'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq','vaoqr',#'vaoqn',
        'vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy','vaoqz','vaqgl','vaqgm','vaqgn']
nina_ensembles = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl','vamrm',
        'vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt','vamru','vamrv','vamrw','vamrx',
        'vamry','vamrz','vaqga','vaqgb','vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']

# Directory of data
directory = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
filename = directory+'vamrk/EHF_heatwaves_ACCESS1-0_vamrk_yearly_summer.nc'

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

# Load control ensebles
for n, ens in enumerate(control_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1-0_'+ens+'_yearly_summer.nc'
    control['hwf'][n], control['hwn'][n], \
    control['hwd'][n], control['hwa'][n], \
    control['hwm'][n], control['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load el nino ensebles
for n, ens in enumerate(nino_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1-0_'+ens+'_yearly_summer.nc'
    elnino['hwf'][n], elnino['hwn'][n], \
    elnino['hwd'][n], elnino['hwa'][n], \
    elnino['hwm'][n], elnino['hwt'][n], _, _ = load_ensemble_hw(filename)
# Load la nina ensebles
for n, ens in enumerate(nina_ensembles):
    filename = directory+ens+'/EHF_heatwaves_ACCESS1-0_'+ens+'_yearly_summer.nc'
    lanina['hwf'][n], lanina['hwn'][n], \
    lanina['hwd'][n], lanina['hwa'][n], \
    lanina['hwm'][n], lanina['hwt'][n], _, _ = load_ensemble_hw(filename)

def plotmap_test(data,name,colours='bwr'):
    fmt = '%1.0f'
    fig = plt.figure()
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
    shade = map_axes.pcolormesh(px,py,data,cmap=colours,vmin=0,vmax=1)
    map_axes.drawparallels(parallels, labels=[True,False,False,False], linewidth=0.5)
    map_axes.drawmeridians(meridians, labels=[False,False,False,True], linewidth=0.5)
    map_axes.drawcoastlines()
    plt.title(name[9:16]+' '+name[0:4]+' '+name[5:8])
    plt.savefig(name, format='eps')
    plt.close()

# Test for normality
msk = control['hwf'].mask[0,...] 
for name, aspect in control.items():
    sksq = np.ones(aspect.shape[-2:])*np.nan
    pval = np.ones(aspect.shape[-2:])*np.nan
    for x in xrange(aspect.shape[2]):
        for y in xrange(aspect.shape[1]):
            if aspect.mask[:,y,x].sum()<8:
                sksq[y,x], pval[y,x] = stats.mstats.normaltest(aspect[:,y,x])
    sksq = np.ma.array(sksq, mask=msk)
    pval = np.ma.array(pval, mask=msk)
    # If the probability (pv) that the sample is normal is less than 0.05 
    # then it gets a value of 0 False
    sig = np.ma.array(pval>=0.05, mask=msk)
    # RdYlGn is reversed so blue is 1 True and red is 0 False.
    plotmap_test(sig, 'normality_'+name+'.eps', colours='bwr_r')

# Test for equal variance
for (cname, caspect), (dname, daspect) in zip(control.items(), elnino.items()):
    ws = np.ones(caspect.shape[-2:])*np.nan
    pval = np.ones(caspect.shape[-2:])*np.nan
    for x in xrange(caspect.shape[2]):
        for y in xrange(caspect.shape[1]):
            if caspect.mask[:,y,x].sum()<8:
                ws[y,x], pval[y,x] = stats.levene(caspect[:,y,x], daspect[:,y,x])
    ws = np.ma.array(ws, mask=msk)
    pval = np.ma.array(pval, mask=msk)
    sig = np.ma.array(pval>=0.05, mask=msk)
    plotmap_test(sig, 'equalvar_'+cname+'_nino.eps', colours='bwr_r')
for (cname, caspect), (dname, daspect) in zip(control.items(), lanina.items()):
    ws = np.ones(caspect.shape[-2:])*np.nan
    pval = np.ones(caspect.shape[-2:])*np.nan
    for x in xrange(caspect.shape[2]):
        for y in xrange(caspect.shape[1]):
            if caspect.mask[:,y,x].sum()<8:
                ws[y,x], pval[y,x] = stats.levene(caspect[:,y,x], daspect[:,y,x])
    ws = np.ma.array(ws, mask=msk)
    pval = np.ma.array(pval, mask=msk)
    sig = np.ma.array(pval>=0.05, mask=msk)
    plotmap_test(sig, 'equalvar_'+cname+'_nina.eps', colours='bwr_r')

# Skewness plots
def plotmap_skew(data,name,colours='bwr'):
    fmt = '%1.0f'
    fig = plt.figure()
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
    shade = map_axes.pcolormesh(px,py,data,cmap=colours,vmin=-2,vmax=2)
    map_axes.drawparallels(parallels, labels=[True,False,False,False], linewidth=0.5)
    map_axes.drawmeridians(meridians, labels=[False,False,False,True], linewidth=0.5)
    map_axes.drawcoastlines()
    cb = plt.colorbar(shade, orientation='horizontal')
    cb.ax.set_xlabel('Skewness')
    plt.title(name[9:16]+' '+name[0:4]+' '+name[5:8])
    plt.tight_layout()
    plt.savefig(name, format='eps')
    plt.close()

for aspect in control:
    plotmap_skew(stats.mstats.skew(control[aspect], axis=0), 'skew_'+aspect+'_control.eps', colours='bwr')
    plotmap_skew(stats.mstats.skew(elnino[aspect], axis=0), 'skew_'+aspect+'_nino.eps', colours='bwr')
    plotmap_skew(stats.mstats.skew(lanina[aspect], axis=0), 'skew_'+aspect+'_nina.eps', colours='bwr')

# Plot PDFs of each aspect for a single gridpoint (Melbourne)
for aspect in control:
    nonan = control[aspect][:,17,5]
    nonan = nonan[nonan.mask==False]
    density = stats.gaussian_kde(nonan)
    xs = np.linspace(0,nonan.max()+5,200)
    plt.plot(xs,density(xs),'k')
    plt.axvline(x=nonan.mean(),color='k')
    plt.axvline(x=np.ma.median(nonan),color='k',linestyle='--')
    nonan = elnino[aspect][:,17,5]
    nonan = nonan[nonan.mask==False]
    density = stats.gaussian_kde(nonan)
    plt.plot(xs,density(xs),'r')
    plt.axvline(x=nonan.mean(),color='r')
    plt.axvline(x=np.ma.median(nonan),color='r',linestyle='--')
    nonan = lanina[aspect][:,17,5]
    nonan = nonan[nonan.mask==False]
    density = stats.gaussian_kde(nonan)
    plt.plot(xs,density(xs),'b')
    plt.axvline(x=nonan.mean(),color='b')
    plt.axvline(x=np.ma.median(nonan),color='b',linestyle='--')
    plt.title('PDF of Melbourne '+aspect)
    plt.savefig('pdf_'+aspect+'.eps', format='eps')
    plt.close()

def plotmaps(data,sig,name,filename,colours='viridis'):
    fmt = '%1.0f'
    if name=='hwf':
        units = 'Days'
        cints = np.arange(-9,10,3)
    elif name=='hwd':
        units = 'Days'
        cints = np.arange(-4,5,1)
    elif name=='hwn':
        units = 'Heatwaves'
        cints = np.arange(-2,2.4,0.4)
        fmt = '%1.1f'
    elif name=='hwa':
        units = '$^{\circ}C^{2}$'
        cints = np.arange(-10,12,2)
    elif name=='hwm':
        units = '$^{\circ}C^{2}$'
        cints = np.arange(-5,6,1)
    elif name=='hwt':
        units = 'Days'
        cints = np.arange(-30,31,10)
    else: 
        units = 'stuff'
    fig = plt.figure()
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
    shade = map_axes.pcolormesh(px,py,data,cmap=colours,vmin=cints[0],vmax=cints[-1])
    mask = map_axes.contourf(x,y,sig, 1, colors='none', hatches=[None,'x'])
    cont = map_axes.contour(x,y,data,levels=cints,colors='k',linewidths=0.5)
    for c in cont.collections:
        if c.get_linestyle() == [(None, None)]:
            continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    plt.clabel(cont, fmt=fmt, fontsize=10)
    cb = plt.colorbar(shade, ticks=cints, orientation='horizontal')
    cb.ax.set_xlabel(units)
    map_axes.drawparallels(parallels, labels=[True,False,False,False], linewidth=0.5)
    map_axes.drawmeridians(meridians, labels=[False,False,False,True], linewidth=0.5)
    map_axes.drawcoastlines()
    plt.title(name)
    plt.tight_layout()
    plt.savefig(filename, format='eps')
    plt.close()


def kruskal_2d(ad,bd):
    hs = np.ones(msk.shape)*np.nan
    hs = np.ma.array(hs,mask=msk)
    pv = hs.copy()
    for y in xrange(ad.shape[1]):
        for x in xrange(bd.shape[2]):
            if hs.mask[y,x]==False:
                hs[y,x], pv[y,x] = stats.kruskal(ad[:,y,x],bd[:,y,x],nan_policy='omit')
    return hs, pv

# Nino diff
for aspect in control:
    # Differene of means for elnino
    ts, pv = stats.ttest_ind(elnino[aspect],control[aspect],axis=0,equal_var=False)
    sig = pv<0.05
    nino_diff = np.ma.mean(elnino[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(nino_diff, sig, aspect, 'mean_'+aspect+'_nino_diff.eps', 'bwr')
    # La nina
    ts, pv = stats.ttest_ind(lanina[aspect],control[aspect],axis=0,equal_var=False)
    sig = pv<0.05
    nina_diff = np.ma.mean(lanina[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(nina_diff, sig, aspect, 'mean_'+aspect+'_nina_diff.eps', 'bwr')
    # el nino vs la nina
    ts, pv = stats.ttest_ind(elnino[aspect],lanina[aspect],axis=0,equal_var=False)
    sig = pv<0.05
    nino_nina_diff = np.ma.mean(elnino[aspect], axis=0)-np.ma.mean(lanina[aspect], axis=0)
    plotmaps(nino_nina_diff, sig, aspect, 'mean_'+aspect+'_nino_nina_diff.eps', 'bwr')

    # Difference of medians nino
    ad = elnino[aspect].data
    bd = control[aspect].data
    hs, pv = kruskal_2d(ad,bd)
    sig = pv<0.05
    nino_diff = np.ma.median(elnino[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(nino_diff, sig, aspect, 'median_'+aspect+'_nino_diff.eps', 'bwr')
    # nina
    ad = lanina[aspect].data
    bd = control[aspect].data
    hs, pv = kruskal_2d(ad,bd)
    sig = pv<0.05
    nina_diff = np.ma.median(lanina[aspect], axis=0)-np.ma.mean(control[aspect], axis=0)
    plotmaps(nina_diff, sig, aspect, 'median_'+aspect+'_nina_diff.eps', 'bwr')
    # nino-nina
    ad = elnino[aspect].data
    bd = lanina[aspect].data
    hs, pv = kruskal_2d(ad,bd)
    sig = pv<0.05
    nino_nina_diff = np.ma.median(elnino[aspect], axis=0)-np.ma.mean(lanina[aspect], axis=0)
    plotmaps(nino_nina_diff, sig, aspect, 'median_'+aspect+'_nino_nina_diff.eps', 'bwr')

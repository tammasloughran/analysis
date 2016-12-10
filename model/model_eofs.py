# -*- coding: utf-8 -*-
import numpy as np
import netCDF4 as nc
from eofs.standard import Eof
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import rotate


def load_heatwaves(filename):
    """load the heatwaves data from a file and select just the australian region.
    
    Argument
    filename -- filename of ensemble yearly heatwave file
    
    Returns
    enshwf, enshwn, enshwd, enshwa, enshwm, enshwt
    """
    hwfile = nc.Dataset(filename, 'r')
    lats = hwfile.variables['lat'][:]
    lons = hwfile.variables['lon'][:]
    enshwf = hwfile.variables['HWF_EHF'][0]
    enshwf = enshwf[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    enshwn = hwfile.variables['HWN_EHF'][0]
    enshwn = enshwn[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    enshwd = hwfile.variables['HWD_EHF'][0]
    enshwd = enshwd[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    enshwa = hwfile.variables['HWA_EHF'][0]
    enshwa = enshwa[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    enshwm = hwfile.variables['HWM_EHF'][0]
    enshwm = enshwm[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    enshwt = hwfile.variables['HWT_EHF'][0]
    enshwt = enshwt[(lats<-10.)&(lats>-44.),:][:,(lons<156.)&(lons>112.)]
    hwfile.close()
    return enshwf, enshwn, enshwd, enshwa, enshwm, enshwt


def nsigpcs(varfrac, north):
    upper = varfrac + north
    lower = varfrac - north
    upper = upper[1:]
    lower = lower[:-1]
    clear = upper<lower
    return np.argmax(clear==False)


# Define Heatwave directories and ensemble codes
hwdir = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
control = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh',
           'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp',
           'vaowq','vaowr','vaows','vaowt','vaowu','vaowv','vaoww','vaowx',
           'vaowy','vaowz','vaqgi','vaqgj','vaqgk','vaowf']
elnino = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg',
          'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq',
          'vaoqr','vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy',
          'vaoqz','vaqgl','vaqgm','vaqgn','vaoqh','vaoqn']
lanina = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl',
          'vamrm','vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt',
          'vamru','vamrv','vamrw','vamrx','vamry','vamrz','vaqga','vaqgb',
          'vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']


# Load netcdf metadata
filename = '/srv/ccrc/data46/z5032520/ehfheatwaves/vamrd/EHF_heatwaves_ACCESS1.3_vamrd_yearly_summer.nc'
example_file = nc.Dataset(filename,'r')
example_data,_,_,_,_,_ = load_heatwaves(filename)
nlatnlon = example_data.shape[-2:]
lats = example_file.variables['lat'][:]
lons = example_file.variables['lon'][:]
lats = lats[(lats<-10.)&(lats>-44.)]
lons = lons[(lons<156.)&(lons>112.)]
example_file.close()


# Allocate heatwave arrays.
nens = len(control) + len(elnino) + len(lanina)
hwf = np.ma.ones((nens,)+nlatnlon)*np.nan
hwn = np.ma.ones((nens,)+nlatnlon)*np.nan
hwd = np.ma.ones((nens,)+nlatnlon)*np.nan
hwa = np.ma.ones((nens,)+nlatnlon)*np.nan
hwm = np.ma.ones((nens,)+nlatnlon)*np.nan
hwt = np.ma.ones((nens,)+nlatnlon)*np.nan


# Load the ensembles and construct a dummy time axis as lanina->control->elnino
i = 0
for ens in lanina:
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    hwf[i,...], hwn[i,...], hwd[i,...], \
    hwa[i,...], hwm[i,...], hwt[i,...] = load_heatwaves(filename)
    i += 1
for ens in control:
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    hwf[i,...], hwn[i,...], hwd[i,...], \
    hwa[i,...], hwm[i,...], hwt[i,...] = load_heatwaves(filename)
    i += 1
for ens in elnino:
    filename = hwdir+ens+'/EHF_heatwaves_ACCESS1.3_'+ens+'_yearly_summer.nc'
    hwf[i,...], hwn[i,...], hwd[i,...], \
    hwa[i,...], hwm[i,...], hwt[i,...] = load_heatwaves(filename)
    i += 1

hwd[(np.logical_not(hwf.mask)&hwd.mask)] = 0
hwa[(np.logical_not(hwf.mask)&hwa.mask)] = 0
hwm[(np.logical_not(hwf.mask)&hwm.mask)] = 0
hwt[(np.logical_not(hwf.mask)&hwt.mask)] = 0

hwasp = [hwf,hwn,hwd,hwa,hwm,hwt]
label = ['HWF','HWN','HWD','HWA','HWM','HWT']
levels = [np.linspace(-10, 10, 11),\
        np.linspace(-2, 2, 11),\
        np.linspace(-5, 5, 11),\
        np.linspace(-20, 20, 11),\
        np.linspace(-8, 8, 11),\
        np.linspace(-50, 50, 11)]

for n,aspect in enumerate(hwasp):
    # Perform the rotated PCA
    coslat = np.cos(np.deg2rad(lats)).clip(0.,1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(aspect, weights=wgts)
    explained_variance = solver.varianceFraction()
    errors = solver.northTest(vfscaled=True)
    #retain = nsigpcs(explained_variance, errors)
    eofs = solver.eofs(eofscaling=2, neofs=4)
    pcs = solver.pcs(pcscaling=1, npcs=4)
    pcs, eofs, R = rotate.do_rotation(pcs, eofs)
    
    # Multiply negative eofs by -1
    for ii in xrange(eofs.shape[0]):
        if eofs[ii].mean()<0.:
            eofs[ii] *= -1
            pcs[...,ii] *= -1

    # Plot some figures
    fig, axes = plt.subplots(nrows=2, ncols=2)
    for i,ax in enumerate(axes.flat):
        if (i>=eofs.shape[0]): continue
        mp = Basemap(ax=ax, projection='cea', llcrnrlon=lons[0], llcrnrlat=lats[0],
                     urcrnrlon=lons[-1], urcrnrlat=lats[-1], resolution='i')
        x,y = np.meshgrid(lons,lats)
        xx,yy = mp(x,y)
        pcol = mp.contourf(xx,yy,eofs[i], levels=levels[n], cmap='seismic')
        cntrs = mp.contour(xx,yy,eofs[i], levels=levels[n], colors='k')
        mp.drawcoastlines()
        mp.drawmeridians([120,130,140,150],linewidth=0.1,labels=[1,0,0,1])
        mp.drawparallels([-15,-20,-25,-30,-35,-40],linewidth=0.1,labels=[0,1,1,0])
    cbar_ax, kw = mpl.colorbar.make_axes([axs for axs in axes.flat], orientation='horizontal')
    cb = plt.colorbar(pcol, cax=cbar_ax, orientation='horizontal')
    plt.title(label[n]+' EOFs')
    plt.savefig(label[n]+'_EOFs.eps',format='eps')
    plt.close()

    for pc in xrange(4):
        if pc>=pcs.shape[-1]: continue
        plt.figure()
        groups = np.zeros((90))
        groups[30:60] = 1
        groups[60:90] = 2
        plt.clf()
        plt.scatter(groups,pcs[...,pc])
        plt.boxplot(pcs[0:30,pc],positions=[0])
        plt.boxplot(pcs[30:60,pc],positions=[1])
        plt.boxplot(pcs[60:90,pc],positions=[2])
        plt.xticks([0,1,2], ['La Nina','Control','El Nino'])
        plt.xlim((-.5,2.5))
        plt.axhline(color='k')
        plt.ylabel(label[n]+' PC'+str(pc+1))
        plt.savefig(label[n]+'_PC'+str(pc+1)+'.png',format='png')
        plt.close()

    # Construct the masks
    masks = np.zeros((4,)+eofs.shape[-2:])
    for msk in xrange(0,4,1):
        sigma = eofs[msk,...].std()
        masks[msk,...] = eofs[msk,...]>sigma*2.0

    # Save the masks to file
    outnc = nc.Dataset(label[n]+'_masks.nc','w')
    outnc.createDimension('lons', masks.shape[2])
    outnc.createDimension('lats', masks.shape[1])
    outnc.createDimension('nmask', masks.shape[0])
    olons = outnc.createVariable('lons', float, dimensions=('lons'))
    olats = outnc.createVariable('lats', float, dimensions=('lats'))
    onmask = outnc.createVariable('nmask', int, dimensions=('nmask'))
    omasks = outnc.createVariable('masks', float, dimensions=('nmask','lats','lons'))
    olons[:] = lons
    olats[:] = lats
    onmask[:] = range(masks.shape[0])
    omasks[:] = masks
    outnc.close()
    np.save(label[n]+'_masks.npy', masks)
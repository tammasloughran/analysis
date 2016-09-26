import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import pdb
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import glob
import pdb


# Ensemble codes
control = ['vamrc','vaowa','vaowb','vaowc','vaowd','vaowe','vaowg','vaowh',
        'vaowi','vaowj','vaowk','vaowl','vaowm','vaown','vaowo','vaowp',
        'vaowq','vaowr','vaows','vaowt','vaowu','vaowv','vaoww','vaowx',
        'vaowy','vaowz','vaqgi','vaqgj','vaqgk']
elnino = ['vamrd','vaoqa','vaoqb','vaoqc','vaoqd','vaoqe','vaoqf','vaoqg',
        'vaoqi','vaoqj','vaoqk','vaoql','vaoqm','vaoqo','vaoqp','vaoqq',
        'vaoqr','vaoqs','vaoqt','vaoqu','vaoqv','vaoqw','vaoqx','vaoqy',
        'vaoqz','vaqgl','vaqgm','vaqgn']
lanina = ['vamre','vamrf','vamrg','vamrh','vamri','vamrj','vamrk','vamrl',
        'vamrm','vamrn','vamro','vamrp','vamrq','vamrr','vamrs','vamrt',
        'vamru','vamrv','vamrw','vamrx','vamry','vamrz','vaqga','vaqgb',
        'vaqgc','vaqgd','vaqge','vaqgf','vaqgg','vaqgh']

# directories
hwdir = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'

# Load Lats and lons
gridnc = nc.Dataset(hwdir+'vamrc/EHF_heatwaves_ACCESS1-0_vamrc_daily.nc','r')
lats = gridnc.variables['lat'][:]
lons = gridnc.variables['lon'][:]
gridnc.close()

# Construct dates
dates = pd.date_range(dt.datetime(2000,1,1),dt.datetime(2001,12,31))


def load_events(filename):
    """Loads heatwave event indicators from an ensemble daily heatwave file.

    Arguments:
    filename

    Returns:
    events
    """
    openfile = nc.Dataset(filename,'r')
    events = openfile.variables['events'][:]
    openfile.close()
    return events


def load_ends(filename):
    """Loads heatwave event indicators from an ensemble daily heatwave file.

    Arguments:
    filename

    Returns:
    events
    """
    openfile = nc.Dataset(filename,'r')
    events = openfile.variables['ends'][:]
    openfile.close()
    return events


def load_pr(filename):
    """Loads pressure data from an ensemble model output file.

    Arguments:
    filename

    Returns:
    pr
    """
    openfile = nc.MFDataset(filename,'r')
    pr = openfile.variables['p_1'][:]
    openfile.close()
    return pr.squeeze()


def plot_pr(data, n, title, filename):
    # Create the map projection
    m = Basemap(projection='mill',
            #boundinglat=-10,lon_0=90,resolution='l')
            llcrnrlon=70.,llcrnrlat=-60.,
            urcrnrlon=230.,urcrnrlat=20.)
    # Grid the lats and lons and feed it to the map instance
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    # Define the range of the pressure levels to contour
    #levs = np.arange(900,1150,5)
    # Plot contours and filled contours
    conts = m.contour(x,y,data/100.,linewidths=0.4,colors='k')
    plt.clabel(conts, fmt='%4.3f', fontsize=10)
    #m.contourf(x,y,data/100.,cmap='viridis')
    # Make the colourbar
    #cbar = plt.colorbar(orientation='horizontal')
    #cbar.set_label('hPa')
    # Plot the significance mask
    #sigmask, p_fdr = bh_fdr.bh_fdr(p, 0.05)
    #sigmask = p<0.05
    #m.contourf(x, y, sigmask, 1, colors='none', hatches=[None,'x'])
    m.drawcoastlines()
    m.drawmeridians(np.arange(lons[0],lons[-1]+10,40.),labels=[1,0,0,1],linewidth=0)
    m.drawparallels(np.arange(-60,1,20.),labels=[1,0,0,1],linewidth=0)
    plt.title(title+' consisting of '+str(n)+'days')
    plt.savefig(filename, format='png')
    #plt.show()
    plt.close()


def select_heatwave_days_pr(ensemble, group):
    """Select the third day of pressure for heatwaves >=5 days duration

    Arguments:
    ensembles -- a strings with the ensemble code
    group -- a string describing the experiment group control elnino or lanina

    Returns:
    composite
    """
    # Load the heatwave events.
    #events = load_events(hwdir+ensemble+'/EHF_heatwaves_ACCESS1-0_'+ensemble+'_daily.nc')
    ends = load_ends(hwdir+ensemble+'/EHF_heatwaves_ACCESS1-0_'+ensemble+'_daily.nc')
    # Exclude events < 5 days
    ends[ends<5] = 0

    # Exlude out of season events; 304=Nov 1st, 455=Apr 1st 
    ends = ends[:455,...]
    ends = ends[304:,...]
    
    # Select the region of interest (this will change at some point)
    # Northeastern Australia
    #reg = ends[:,:,75:80]
    #reg = nea[:,51:56,:]
    # Southeastern Australia
    reg = ends[:,:,75:80]
    reg = reg[:,42:46,:]
    
    # Construct a time index indicating if heatwaves occur within this region.
    index = np.max(np.max(reg,axis=2), axis=1)
    for i in np.where(index>0)[0]:
        try:
            for ii in xrange(0,int(index[i]-1)):
                index[i+ii] += 1
        except IndexError:
            continue
    index = index>0

    # Select only the third day of heatwave
    #index = reg.sum(axis=2)
    #index = index.sum(axis=1)
    #index = np.diff(index) # shape - 1
    #index[index<0] = 0
    #index = np.append([0],index) # restore original shape
    #for i in np.where(index>0)[0]:
    #    index[i] = 0
    #    index[i+2] = 1

    # Load the pressure data
    prfiles = [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2000-11.nc']
    prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2000-12.nc']
    prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-01.nc']
    prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-02.nc']
    prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-03.nc']
    pr = load_pr(prfiles)
    
    # Select only the third heatwave days using index
    pr = pr[index,...]
    
    # Calculate anomalies from climatology
    return pr - pclim


def construct_composite(ensembles, group):
    """Construct a composite of all heatwaves within the northeast australian 
    region for all ensembles in a group.

    Arguments:
    ensembles -- a list of strings with the ensemble codes
    group -- a string describing the experiment group control elnino or lanina

    Returns:
    composite
    """
    composite = np.ones((len(ensembles), len(lats), len(lons)))
    for i, ens in enumerate(ensembles):
        # Load the heatwave events.
        events = load_events(hwdir+ens+'/EHF_heatwaves_ACCESS1-0_'+ens+'_daily.nc')
    
        # Select the northeastern Australian region
        nea = events[:,:,75:80]
        nea = nea[:,51:56,:]

        # Construct a time index indicating if heatwaves occur within this region.
        index = nea.sum(axis=2)
        index = index.sum(axis=1)
        index = index>0

        # Load the pressure data
        pr = load_pr(modeldir+group+'/'+ens+'/'+ens+'a.pe????-??.nc')
        pr = np.delete(pr, [59], axis=0) # Remove feb 29th. It doesnt exist in events

        # Select only the heatwave days using index
        pr = pr[index,...]/100.

        # Composite of heatwave days
        composite[i,...] = np.mean(pr, axis=0)
    composite = np.mean(composite, axis=0)
    return composite


def plotmaps(data,name,filename,colours='viridis'):
    fmt = '%1.0f'
    parallels = np.arange(-90., 120., 10.)
    #meridians = np.arange(0., 360., 10.)
    map_axes = Basemap(projection='robin',lon_0=180.,resolution='l')
    map_axes.drawparallels(parallels, linewidth=0.3, \
                    labels=[True,False,False,False], fontsize=8, xoffset=0.6)
    xx, yy = np.meshgrid(lons, lats)
    x, y = map_axes(xx,yy)
    xx = xx - 0.5 # The data projection is slightly off compared to the coastlines.
    yy = yy - 0.5
    px, py = map_axes(xx,yy)
    #data = np.ma.array(data,mask=np.isnan(data))
    shade = map_axes.pcolormesh(px,py,data,cmap=colours)
    #shade = map_axes.pcolormesh(px,py,data,cmap=colours,vmin=cints[0],vmax=cints[-1])
    #mask = map_axes.contourf(x,y,sig, 1, colors='none', hatches=[None,'xx'])
    cont = map_axes.contour(x,y,data,colors='k',linewidths=0.5)
    #cont = map_axes.contour(x,y,data,levels=cints,colors='k',linewidths=0.5)
    for c in cont.collections:
        if c.get_linestyle() == [(None, None)]:
            continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    plt.clabel(cont, fmt=fmt, fontsize=10)
    cb = map_axes.colorbar(shade, location='bottom', pad=0.25)
    #cb.ax.set_xlabel(units)
    #map_axes.drawparallels(parallels, labels=[True,False,False,False], linewidth=0.5)
    #map_axes.drawmeridians(meridians, labels=[False,False,False,True], linewidth=0.5)
    map_axes.drawcoastlines()
    plt.title(name)
    plt.tight_layout()
    #plt.savefig(filename, format='png')
    plt.show()
    plt.close()


if __name__=='__main__':
    ncs = nc.MFDataset('/srv/ccrc/data46/z5032520/modelout/ACCESS/vamrb/vamrba.pa19*')
    pclim = ncs.variables['p'][-12*30:]
    ncs.close()
    pclim = np.squeeze(np.mean(pclim, axis=0))
    comp = np.zeros(pclim.shape)
    n = 0
    for ens in control:
        print ens
        # Select heatwave days
        pr_days = select_heatwave_days_pr(ens,'control')
        
        # plot the maps
        #for i in xrange(pr_days.shape[0]):
        #    plot_pr(pr_days[i,...], '%s %s'%(ens, i), '%s_%s.eps'%(ens, i))
        if pr_days.shape[0]>0:
            #plot_pr(pr_days.mean(axis=0)-pclim, '%s composite'%(ens), '%s_composite.png'%(ens))
            n += 1
            ndays = n + pr_days.shape[0]
            comp = comp + np.mean(pr_days, axis=0)
    ccomp = comp/float(n)
    plot_pr(ccomp, ndays, 'control composite', 'control_composite.png')
    comp = np.zeros(pclim.shape)
    n = 0
    for ens in elnino:
        print ens
        # Select heatwave days
        pr_days = select_heatwave_days_pr(ens,'elnino')

        # plot the maps
        #for i in xrange(pr_days.shape[0]):
        #    plot_pr(pr_days[i,...], '%s %s'%(ens, i), '%s_%s.eps'%(ens, i))
        if pr_days.shape[0]>0:
            #plot_pr(pr_days.mean(axis=0)-pclim, '%s composite'%(ens), '%s_composite.png'%(ens))
            n += 1
            ndays = n + pr_days.shape[0]
            comp = comp + np.mean(pr_days, axis=0)
    ocomp = comp/float(n)
    plot_pr(ocomp, ndays, 'El Nino composite', 'elnino_composite.png')
    comp = np.zeros(pclim.shape)
    n = 0
    for ens in lanina:
        print ens
        # Select heatwave days
        pr_days = select_heatwave_days_pr(ens,'lanina')

        # plot the maps
        #for i in xrange(pr_days.shape[0]):
        #    plot_pr(pr_days[i,...], '%s %s'%(ens, i), '%s_%s.eps'%(ens, i))
        if pr_days.shape[0]>0:
            #plot_pr(pr_days.mean(axis=0)-pclim, '%s composite'%(ens), '%s_composite.png'%(ens))
            n += 1  
            ndays = n + pr_days.shape[0]
            comp = comp + np.mean(pr_days, axis=0)
    acomp = comp/float(n)
    plot_pr(acomp, ndays, 'La Nina composite', 'lanina_composite.png')

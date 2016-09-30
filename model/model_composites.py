import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pdb


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


def load_time(filename):
    """Loads a dates from a netcdf file.
    """
    openfile = nc.Dataset(filename,'r')
    times = openfile.variables['time']
    dates = nc.num2date(times[:],units=times.units,calendar=times.calendar)
    drange = pd.period_range(str(dates[0]),str(dates[-1]),freq='D')
    if times.calendar=='365_day':
        drange = drange[(drange.month!=2)|(drange.day!=29)]
    return drange


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
    '''Plots a map pf pressure for Australia and the surrounding oceans
    
    Arguments:
    data -- tThe 2D pressure field to plot
    n -- the number of ensembles that contributed to the composite
    title -- Plot title
    filename -- the output filename
    '''
    # Create the map projection
    m = Basemap(projection='mill',
            #boundinglat=-10,lon_0=90,resolution='l')
            llcrnrlon=70.,llcrnrlat=-60.,
            urcrnrlon=230.,urcrnrlat=20.)
    # Grid the lats and lons and feed it to the map instance
    lns,lts = np.meshgrid(lons,lats)
    x,y = m(lns,lts)
    # Plot contours
    conts = m.contour(x,y,data/100.,linewidths=0.4,colors='k')
    # Plot labels on contours 4 spaces with 2 decimal points
    plt.clabel(conts, fmt='%5.1f', fontsize=10)
    # Map details
    m.drawcoastlines()
    m.drawmeridians(np.arange(lons[0],lons[-1]+10,40.),labels=[1,0,0,1],linewidth=0)
    m.drawparallels(np.arange(-60,1,20.),labels=[1,0,0,1],linewidth=0)
    plt.title(title+' consisting of '+str(n)+' days')
    plt.savefig(filename, format='png')
    plt.close()


def select_heatwave_days_pr(ensemble, group, hwdir, modeldir):
    """Select the third day of pressure for heatwaves >=5 days duration

    Arguments:
    ensembles -- a strings with the ensemble code
    group -- a string describing the experiment group control elnino or lanina

    Returns:
    composite
    """
    # Load the heatwave events.
    ends = load_ends(hwdir+ensemble+'/EHF_heatwaves_ACCESS1.3_'+ensemble+'_daily.nc')
    
    # Exclude events < 5 days
    ends[ends<5] = 0

    # Exlude out of season events.
    dates = load_time(hwdir+ensemble+'/EHF_heatwaves_ACCESS1.3_'+ensemble+'_daily.nc')
    nov1 = np.where((dates.month==11)&(dates.day==1))[0][0]
    mar31 = np.where((dates.month==3)&(dates.day==31))[0][1]
    ends = ends[:mar31+1,...]
    ends = ends[nov1:,...]
    
    # Select the region of interest (this will change at some point)
    # Northwestern Australia
    reg = ends[:,:,67:73]
    reg = reg[:,56:61,:]
    # Northeastern Australia
    #reg = ends[:,:,75:80]
    #reg = reg[:,51:56,:]
    # Southeastern Australia
    #reg = ends[:,:,75:80]
    #reg = reg[:,42:46,:]
    
    # Construct a time index indicating if heatwaves occur within this region.
    index = np.max(np.max(reg,axis=2), axis=1)
    for i in np.where(index>0)[0]:
        try:
            for ii in xrange(0,int(index[i]-1)):
                index[i+ii] += 1
        except IndexError:
            continue
    index = index>0

    # Load the pressure data. Nov-Mar and subtract climatology
    pr = np.empty((0,)+pclim.shape[-2:])
    for i,month in enumerate(['11','12','01','02','03']):
        if month=='11' or month=='12': 
            year = '2000'
        else:
            year = '2001'
        prfile = modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe'+year+'-'+month+'.nc'
        prdata = load_pr(prfile)
        climp = pclim[int(month)-1,...]
        pr2 = prdata - climp
        pr = np.concatenate((pr, pr2), axis=0)

    #prfiles = [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2000-11.nc']
    #prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2000-12.nc']
    #prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-01.nc']
    #prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-02.nc']
    #prfiles += [modeldir+group+'/'+ensemble+'/'+ensemble+'a.pe2001-03.nc']
    #pr = load_pr(prfiles)
    
    # Select only the third heatwave days using index
    pr = pr[index,...]
    
    # Calculate anomalies from climatology
    return pr #- pclim


def construct_composite(ensembles, group, hwdir, modeldir):
    """Construct a composite of all heatwaves within a region for all 
    ensembles in a group.

    Arguments:
    ensembles -- a list of strings with the ensemble codes
    group -- a string describing the experiment group control elnino or lanina

    Returns:
    composite -- the composited field
    ndays -- the number of heatwave days that constitute the composite
    n -- the number of ensembles with heatwave days.
    """
    comp = np.zeros(pclim.shape[-2:])
    n = 0
    ndays = 0
    for ens in ensembles:
        print ens
        # Select heatwave days
        pr_days = select_heatwave_days_pr(ens,group,hwdir,modeldir)
        # If the ensemble has heatwaves then put it in the composite
        if pr_days.shape[0]>0:
            n += 1
            ndays += pr_days.shape[0]
            comp += np.mean(pr_days, axis=0)
    comp = comp/float(n) # Each ensemble has equal weighting
    return comp, ndays, n


if __name__=='__main__':
    # Ensemble codes
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

    # directories
    hwdir = '/srv/ccrc/data46/z5032520/ehfheatwaves/'
    modeldir = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'

    # Load Lats and lons
    gridnc = nc.Dataset(hwdir+'vamrc/EHF_heatwaves_ACCESS1.3_vamrc_daily.nc','r')
    lats = gridnc.variables['lat'][:]
    lons = gridnc.variables['lon'][:]
    gridnc.close()

    # Construct dates
    dates = pd.date_range(dt.datetime(2000,1,1),dt.datetime(2001,12,31))   

    # Load climatology    
    ncs = nc.MFDataset('/srv/ccrc/data46/z5032520/modelout/ACCESS/vamrb/vamrba.pa19*')
    prc = np.squeeze(ncs.variables['p'][-12*30:])
    ncs.close()
    date_range = pd.period_range(start='1970-01-01',end='1999-12-31',freq='M')
    pclim = np.zeros((12,)+prc.shape[-2:])
    for month in xrange(1,13,1):
        pclim[month-1,...] = prc[date_range.month==month].mean(axis=0)
    
    # Control composite
    ccomp, ndays, n = construct_composite(control,'control', hwdir, modeldir)
    plot_pr(ccomp, ndays, 'control composite', 'control_composite.png')
    
    # El Nino composite
    ocomp, ndays, n = construct_composite(elnino,'elnino', hwdir, modeldir)
    plot_pr(ocomp, ndays, 'El Nino composite', 'elnino_composite.png')
    
    # La Nina composite
    acomp, ndays, n = construct_composite(lanina,'lanina', hwdir, modeldir)
    plot_pr(acomp, ndays, 'La Nina composite', 'lanina_composite.png')
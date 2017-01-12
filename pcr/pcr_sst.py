import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt


# Summer sstaessure anomalies from HadISST
sstfile = ('/media/Jupiter/observations/HadISST/sst/indo-pacific-so_SST.nc')
hadisst = Dataset(sstfile,'r')
sst = hadisst.variables['sst'][:]
sst_lats = hadisst.variables['lat'][:]
sst_lons = hadisst.variables['lon'][:]
sst_time = hadisst.variables['time'][:]
sst_dates = np.array([dt.datetime(1,1,1) + dt.timedelta(hours=hrs) 
        for hrs in sst_time])
sst_years = np.array([sst_dates[i].year 
        for i in np.arange(sst_dates.size)])
sst_months = np.array([sst_dates[i].month 
        for i in np.arange(sst_dates.size)])
# Slice the desired period from 1910-2013
period = np.where(sst_years>=1911)[0]
sst_dates = sst_dates[period]
sst_months = sst_months[period]
sst_years = sst_years[period]
sst_period = sst[period,...]

# Select the desired months to average over
sel = (sst_months==1)|(sst_months==2)|(sst_months==3)|(sst_months==11)|(sst_months==12)
sst_sel = sst_period[sel,...]

# Calculate annual summer means
sst_annual = np.ma.zeros((102,)+sst_sel.shape[-2:])
for year in range(102):
    sst_annual[year,...] = sst_sel[3+5*year:3+5+5*year,...].mean(axis=0)
ssta = sst_annual
# Calculate anomalies with respect to 1961-1990 base period
#years = np.arange(1911,2013)
#base = np.where((years>=1961)&(years<=1990))[0]
#ave = sst_annual[base,...].mean(axis=0)
#ssta = sst_annual - ave

direct = '/home/nfs/z5032520/analysis/eofs/hw_metrics/'
# Principle components from PCA of heatwaves
hwn_pcs = np.load(direct+'HWN_rotated_pcs.npy')[:]
hwf_pcs = np.load(direct+'HWF_rotated_pcs.npy')[:]
hwd_pcs = np.load(direct+'HWD_rotated_pcs.npy')[:]
hwa_pcs = np.load(direct+'HWA_rotated_pcs.npy')[:]
hwm_pcs = np.load(direct+'HWM_rotated_pcs.npy')[:]
hwt_pcs = np.load(direct+'HWT_rotated_pcs.npy')[:]

# Make a 2D regression function
def linregress_2D(x,y):
    """linregress_2D does the same thing as scipy.stats.linregress but accepts
    1D independant variable x and a 3D dependant variable y. Regression is 
    performed along the first axis.
    """
    yshape = y.shape
    space = yshape[1]*yshape[2]
    y = y.reshape(yshape[0], space)
    slope = np.ones(space)*np.nan
    intercept = np.ones(space)*np.nan
    correlation = np.ones(space)*np.nan
    p = np.ones(space)*np.nan
    error = np.ones(space)*np.nan
    for i in xrange(0, space):
        if y.mask[0,i]: continue
        slope[i], intercept[i], correlation[i], p[i], error[i] = \
                linregress(x,y[:,i])
    slope = slope.reshape(yshape[1:])
    intercept = intercept.reshape(yshape[1:])
    correlation = correlation.reshape(yshape[1:])
    p = p.reshape(yshape[1:])
    error = error.reshape(yshape[1:])
    return slope, intercept, correlation, p, error

def plot_pcr(slope, p, title, filename):
    import math
    import bh_fdr
    plt.figure()
    # Create the map projection
    m = Basemap(projection='mill',
            llcrnrlon=sst_lons[0],llcrnrlat=sst_lats[-1],
            urcrnrlon=sst_lons[-1],urcrnrlat=sst_lats[0])
    # Grid the lats and lons and feed it to the map instance
    lns,lts = np.meshgrid(sst_lons,sst_lats)
    x,y = m(lns,lts)
    # Define the range of the sstaessure levels to contour
    levs = np.arange(-0.6,0.65,0.1)
    # Plot contours and filled contours
    m.contour(x,y,slope,linewidths=0.4,colors='k',levels=levs)
    fcont = m.contourf(x,y,slope,cmap='bwr',levels=levs)
    # Make the colourbar
    cbar = plt.colorbar(fcont, orientation='horizontal')
    cbar.set_label('$^\circ$C')
    # Plot the significance mask
    sigmask, p_fdr = bh_fdr.bh_fdr(p, 0.05)
    #sigmask = p<0.05
    m.contourf(x, y, sigmask, 1, colors='none', hatches=[None,'x'])
    m.drawcoastlines()
    m.drawmeridians(np.arange(sst_lons[0],sst_lons[-1]+10,50.),labels=[1,0,0,1],linewidth=0)
    m.drawparallels(np.arange(sst_lats[-1]+10-0.5,sst_lats[0],20.),labels=[1,0,0,1],linewidth=0)
    plt.title(title)
    plt.savefig(filename, format='eps')
    plt.close()

def plot_pcr_add(slope, p, title, filename, axes):
    import math
    import bh_fdr
    # Create the map projection
    m = Basemap(ax=axes, projection='mill',
            llcrnrlon=sst_lons[0],llcrnrlat=sst_lats[-1],
            urcrnrlon=sst_lons[-1],urcrnrlat=sst_lats[0])
    # Grid the lats and lons and feed it to the map instance
    lns,lts = np.meshgrid(sst_lons,sst_lats)
    x,y = m(lns,lts)
    # Define the range of the sstaessure levels to contour
    levs = np.arange(-0.6,0.65,0.1)
    # Plot contours and filled contours
    m.contour(x,y,slope,linewidths=0.4,colors='k',levels=levs)
    fcont = m.contourf(x,y,slope,cmap='bwr',levels=levs)
    #cbar_ax, kw = mpl.colorbar.make_axes(axes, orientation='horizontal')
    cb = m.colorbar(fcont, pad=0.3, location='bottom')
    cb.set_label('$^\circ$C')
    # Plot the significance mask
    sigmask, p_fdr = bh_fdr.bh_fdr(p, 0.05)
    #sigmask = p<0.05
    m.contourf(x, y, sigmask, 1, colors='none', hatches=[None,'x'])
    m.drawcoastlines()
    m.drawmeridians(np.arange(sst_lons[0],sst_lons[-1]+10,50.),labels=[1,0,0,1],linewidth=0)
    m.drawparallels(np.arange(sst_lats[-1]+10-0.5,sst_lats[0],20.),labels=[1,0,0,1],linewidth=0)
    plt.title(title)
    #plt.savefig(filename, format='eps')
    #plt.close()


# Perform regression and plot
import matplotlib.pyplot as plt
final_fig, final_axes = plt.subplots(nrows=3,ncols=1,figsize=(12,15))
#PC1
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,0], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWN PC1', 'regres_ssta_HWNpc1')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,0], ssta)
plot_pcr(slope, p, 'a) Linear regression of summer SST over HWF PC1', 'regres_ssta_HWFpc1')
first=True
plot_pcr_add(slope, p, 'a) Linear regression of summer SST over HWF PC1', 'regres_ssta_HWFpc1', final_axes[0])
first=False
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,0], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWD PC1', 'regres_ssta_HWDpc1')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,0], ssta)
plot_pcr(slope, p, 'c) Linear regression of summer SST over HWA PC1', 'regres_ssta_HWApc1')
plot_pcr_add(slope, p, 'c) Linear regression of summer SST over HWA PC1', 'regres_ssta_HWApc1', final_axes[2])
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,0], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWM PC1', 'regres_ssta_HWMpc1')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,0], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWT PC1', 'regres_ssta_HWTpc1')
#PC2
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,1], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWN PC2', 'regres_ssta_HWNpc2')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,1], ssta)
plot_pcr(slope, p, 'b) Linear regression of summer SST over HWF PC2', 'regres_ssta_HWFpc2')
plot_pcr_add(slope, p, 'b) Linear regression of summer SST over HWF PC2', 'regres_ssta_HWFpc2', final_axes[1])
final_fig.savefig('final_sst.eps',format='eps')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,1], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWD PC2', 'regres_ssta_HWDpc2')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,1], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWA PC2', 'regres_ssta_HWApc2')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,1], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWM PC2', 'regres_ssta_HWMpc2')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,1], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWT PC2', 'regres_ssta_HWTpc2')
#PC3
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,2], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWN PC3', 'regres_ssta_HWNpc3')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,2], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWF PC3', 'regres_ssta_HWFpc3')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,2], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWD PC3', 'regres_ssta_HWDpc3')
#slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,2], ssta)
#plot_pcr(slope, p, 'Linear regression of summer SSTA over HWA rPC3', 'regres_ssta_HWApc3')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,2], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWM PC3', 'regres_ssta_HWMpc3')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,2], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWT PC3', 'regres_ssta_HWTpc3')
#PC4
#slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,3], ssta)
#plot_pcr(slope, p, 'Linear regression of summer SSTA over HWN rPC4', 'regres_ssta_HWNpc4')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,3], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWF PC4', 'regres_ssta_HWFpc4')
#slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,3], ssta)
#plot_pcr(slope, p, 'Linear regression of summer SSTA over HWD rPC4', 'regres_ssta_HWDpc4')
#slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,3], ssta)
#plot_pcr(slope, p, 'Linear regression of summer SSTA over HWA rPC4', 'regres_ssta_HWApc4')
#slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,3], ssta)
#plot_pcr(slope, p, 'Linear regression of summer SSTA over HWM rPC4', 'regres_ssta_HWMpc4')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,3], ssta)
plot_pcr(slope, p, 'Linear regression of summer SST over HWT PC4', 'regres_ssta_HWTpc4')

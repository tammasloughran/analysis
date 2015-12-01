import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# Summer pressure anomalies from 20CR
prfile = ('/home/nfs/z5032520/analysis/dataproc/summer_mslp_1911-2011.nc')
prnc = Dataset(prfile, 'r')
pr = prnc.variables['mslpa'][:]
pr_lats = prnc.variables['lat'][:]
pr_lons = prnc.variables['lon'][:]
ulat = np.where(pr_lats>-60.)[0].max()
llat = np.where(pr_lats<30.)[0].min()
llon = np.where(pr_lons>40.)[0].min()
ulon = np.where(pr_lons<260.)[0].max()
pr = pr[:,llat:ulat,llon:ulon]
pr_lons = pr_lons[llon:ulon]
pr_lats = pr_lats[llat:ulat]

# Principle components from PCA of heatwaves
hwn_pcs = np.load('HWN_rotated_pcs.npy')[:-1]
hwf_pcs = np.load('HWF_rotated_pcs.npy')[:-1]
hwd_pcs = np.load('HWD_rotated_pcs.npy')[:-1]
hwa_pcs = np.load('HWA_rotated_pcs.npy')[:-1]
hwm_pcs = np.load('HWM_rotated_pcs.npy')[:-1]
hwt_pcs = np.load('HWT_rotated_pcs.npy')[:-1]

# Make a 2D regression function
def linregress_2D(x,y):
    """linregress_2D does the same thing as scipy.stats.linregress but accepts
    1D independant variable x and a 3D dependant variable y. Regression is 
    performed along the first axis.
    """
    yshape = y.shape
    space = yshape[1]*yshape[2]
    y = y.reshape(yshape[0], space)
    slope = np.empty(space)
    intercept = np.empty(space)
    correlation = np.empty(space)
    p = np.empty(space)
    error = np.empty(space)
    for i in xrange(0, space):
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
    # Create the map projection
    m = Basemap(projection='mill',
            llcrnrlon=pr_lons[0],llcrnrlat=pr_lats[-1],
            urcrnrlon=pr_lons[-1],urcrnrlat=pr_lats[0])
    # Grid the lats and lons and feed it to the map instance
    lns,lts = np.meshgrid(pr_lons,pr_lats)
    x,y = m(lns,lts)
    # Define the range of the pressure levels to contour
    levs = np.arange(-120,130,20)
    # Plot contours and filled contours
    m.contour(x,y,slope,linewidths=0.4,colors='k',levels=levs)
    m.contourf(x,y,slope,cmap='bwr',levels=levs)
    # Make the colourbar
    cbar = plt.colorbar(orientation='horizontal')
    cbar.set_label('Pa')
    # Plot the significance mask
    sigmask, p_fdr = bh_fdr.bh_fdr(p, 0.05)
    #sigmask = p<0.05
    m.contourf(x, y, sigmask, 1, colors='none',hatches=[None,'.'])
    m.drawcoastlines()
    m.drawmeridians(np.arange(pr_lons[0],pr_lons[-1]+10,40.),labels=[1,0,0,1],linewidth=0,fontsize=10)
    m.drawparallels(np.arange(pr_lats[-1],pr_lats[0],30.),labels=[1,0,0,1],linewidth=0,fontsize=10)
    plt.title(title)
    plt.savefig(filename ,format='eps')
    plt.close()

# Perform regression and plot
#PC1
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN rPC1', 'regres_pra_HWNpc1.eps')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF rPC1', 'regres_pra_HWFpc1.eps')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD rPC1', 'regres_pra_HWDpc1.eps')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC1', 'regres_pra_HWApc1.eps')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM rPC1', 'regres_pra_HWMpc1.eps')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT rPC1', 'regres_pra_HWTpc1.eps')
#PC2
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN rPC2', 'regres_pra_HWNpc2.eps')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF rPC2', 'regres_pra_HWFpc2.eps')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD rPC2', 'regres_pra_HWDpc2.eps')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC2', 'regres_pra_HWApc2.eps')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM rPC2', 'regres_pra_HWMpc2.eps')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT rPC2', 'regres_pra_HWTpc2.eps')
#PC3
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN rPC3', 'regres_pra_HWNpc3.eps')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF rPC3', 'regres_pra_HWFpc3.eps')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD rPC3', 'regres_pra_HWDpc3.eps')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC3', 'regres_pra_HWApc3.eps')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM rPC3', 'regres_pra_HWMpc3.eps')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT rPC3', 'regres_pra_HWTpc3.eps')
#PC4
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN rPC4', 'regres_pra_HWNpc4.eps')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF rPC4', 'regres_pra_HWFpc4.eps')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD rPC4', 'regres_pra_HWDpc4.eps')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC4', 'regres_pra_HWApc4.eps')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM rPC4', 'regres_pra_HWMpc4.eps')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT rPC4', 'regres_pra_HWTpc4.eps')

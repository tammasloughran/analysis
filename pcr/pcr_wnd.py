import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# Summer pressure anomalies from 20CR
prfile = ('/home/nfs/z5032520/analysis/dataproc/summer_mslp_1911-2011.nc')
prnc = Dataset(prfile, 'r')
pr = prnc.variables['mslp'][:]
pr_lats = prnc.variables['lat'][:]
pr_lons = prnc.variables['lon'][:]
ulat = np.where(pr_lats>-70.)[0].max()
llat = np.where(pr_lats<30.)[0].min()
llon = np.where(pr_lons>40.)[0].min()
ulon = np.where(pr_lons<260.)[0].max()
pr = np.flipud(pr[:,llat:ulat,llon:ulon])
pr_lons = pr_lons[llon:ulon]
pr_lats = np.flipud(pr_lats[llat:ulat])


# Summer wind anomalies from 20CR
ufile = ('/home/nfs/z5032520/analysis/dataproc/summer_uwnd_1911-2011.nc')
vfile = ('/home/nfs/z5032520/analysis/dataproc/summer_vwnd_1911-2011.nc')
unc = Dataset(ufile, 'r')
vnc = Dataset(vfile, 'r')
u = unc.variables['uwnd'][:]
v = vnc.variables['vwnd'][:]
u_lats = unc.variables['lat'][:]
u_lons = unc.variables['lon'][:]
ulat = np.where(u_lats>-70.)[0].max()
llat = np.where(u_lats<30.)[0].min()
llon = np.where(u_lons>40.)[0].min()
ulon = np.where(u_lons<260.)[0].max()
u = np.flipud(u[:,llat:ulat,llon:ulon])
v = np.flipud(v[:,llat:ulat,llon:ulon])
u_lons = u_lons[llon:ulon]
u_lats = np.flipud(u_lats[llat:ulat])


direct = '/home/nfs/z5032520/analysis/eofs/hw_metrics/'
# Principle components from PCA of heatwaves
hwn_pcs = np.load(direct+'HWN_rotated_pcs.npy')[:-1]
hwf_pcs = np.load(direct+'HWF_rotated_pcs.npy')[:-1]
hwd_pcs = np.load(direct+'HWD_rotated_pcs.npy')[:-1]
hwa_pcs = np.load(direct+'HWA_rotated_pcs.npy')[:-1]
hwm_pcs = np.load(direct+'HWM_rotated_pcs.npy')[:-1]
hwt_pcs = np.load(direct+'HWT_rotated_pcs.npy')[:-1]

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
            #boundinglat=-10,lon_0=90,resolution='l')
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
    #uproj,vproj,xx,yy = \
    #        m.transform_vector(slope2,slope3,u_lons,u_lats,slope2.shape[1],slope2.shape[0],returnxy=True,masked=True)
    #Q = m.quiver(xx,yy,uproj,vproj,scale=700)
    #qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
    # Make the colourbar
    cbar = plt.colorbar(orientation='horizontal')
    cbar.set_label('Pa')
    # Plot the significance mask
    sigmask, p_fdr = bh_fdr.bh_fdr(p, 0.05)
    #sigmask = p<0.05
    m.contourf(x, y, sigmask, 1, colors='none', hatches=[None,'x'])
    m.drawcoastlines()
    m.drawmeridians(np.arange(pr_lons[0],pr_lons[-1]+10,40.),labels=[1,0,0,1],linewidth=0)
    m.drawparallels(np.arange(-60,1,20.),labels=[1,0,0,1],linewidth=0)
    plt.title(title)
    plt.savefig(filename, format='eps')
    plt.close()

# Perform regression and plot
#PC1
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,0], pr)
#slope2, intercept, correlation, _, error = linregress_2D(hwn_pcs[:,0], u)
#slope3, intercept, correlation, _, error = linregress_2D(hwn_pcs[:,0], v)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN PC1', 'regres_pr_HWNpc1')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,0], pr)
plot_pcr(slope, p, 'a) Linear regression of summer MSLP over HWF PC1', 'regres_pr_HWFpc1')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD PC1', 'regres_pr_HWDpc1')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,0], pr)
plot_pcr(slope, p, 'c) Linear regression of summer MSLP over HWA PC1', 'regres_pr_HWApc1')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM PC1', 'regres_pr_HWMpc1')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,0], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT PC1', 'regres_pr_HWTpc1')
#PC2
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN PC2', 'regres_pr_HWNpc2')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,1], pr)
plot_pcr(slope, p, 'b) Linear regression of summer MSLP over HWF PC2', 'regres_pr_HWFpc2')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD PC2', 'regres_pr_HWDpc2')
slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA PC2', 'regres_pr_HWApc2')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM PC2', 'regres_pr_HWMpc2')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,1], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT PC2', 'regres_pr_HWTpc2')
#PC3
slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN PC3', 'regres_pr_HWNpc3')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF PC3', 'regres_pr_HWFpc3')
slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD PC3', 'regres_pr_HWDpc3')
#slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,2], pr)
#plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC3', 'regres_pr_HWApc3')
slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM PC3', 'regres_pr_HWMpc3')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,2], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT PC3', 'regres_pr_HWTpc3')
#PC4
#slope, intercept, correlation, p, error = linregress_2D(hwn_pcs[:,3], pr)
#plot_pcr(slope, p, 'Linear regression of summer MSLP over HWN rPC4', 'regres_pr_HWNpc4')
slope, intercept, correlation, p, error = linregress_2D(hwf_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWF PC4', 'regres_pr_HWFpc4')
#slope, intercept, correlation, p, error = linregress_2D(hwd_pcs[:,3], pr)
#plot_pcr(slope, p, 'Linear regression of summer MSLP over HWD rPC4', 'regres_pr_HWDpc4')
#slope, intercept, correlation, p, error = linregress_2D(hwa_pcs[:,3], pr)
#plot_pcr(slope, p, 'Linear regression of summer MSLP over HWA rPC4', 'regres_pr_HWApc4')
#slope, intercept, correlation, p, error = linregress_2D(hwm_pcs[:,3], pr)
#plot_pcr(slope, p, 'Linear regression of summer MSLP over HWM rPC4', 'regres_pr_HWMpc4')
slope, intercept, correlation, p, error = linregress_2D(hwt_pcs[:,3], pr)
plot_pcr(slope, p, 'Linear regression of summer MSLP over HWT PC4', 'regres_pr_HWTpc4')

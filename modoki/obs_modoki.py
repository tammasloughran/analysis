import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
#plt.rc('text', usetex=False)
plt.rc('font', family='sans')
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats
import matplotlib

# Load data
#hwfile = nc.Dataset('/home/nfs/z5032520/analysis/modoki/EHF_heatwaves_20CRV2_1901-2012_yearly_summer.nc', 'r')
#hwfile = nc.Dataset('/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/EHF_heatwaves____yearly_summer.nc','r')
hwfile = nc.Dataset('EHF_heatwaves_AWAP_bp1961-2010_yearly_summer.nc','r')
hwf = hwfile.variables['HWF_EHF'][:]
hwf.data[hwf.mask==True] = np.nan
hwd = hwfile.variables['HWD_EHF'][:]
hwd.data[hwd.mask==True] = np.nan
hwn = hwfile.variables['HWN_EHF'][:]
hwn.data[hwn.mask==True] = np.nan
hwm = hwfile.variables['HWM_EHF'][:]
hwm.data[hwm.mask==True] = np.nan
hwa = hwfile.variables['HWA_EHF'][:]
hwa.data[hwa.mask==True] = np.nan
hwt = hwfile.variables['HWT_EHF'][:]
hwt.data[hwt.mask==True] = np.nan
time = hwfile.variables['time'][:]
lats = hwfile.variables['lat'][:]
lons = hwfile.variables['lon'][:]

base = (time>=1961)&(time<=2010)

# Define modoki years
#modoki_years = [1923, 1929, 1940, 1946, 1958, 1963, 1977, 1986,
#                1990, 1991, 1992, 1994, 2002, 2004, 2009]
modoki_years = [1986, 1990, 1991, 1992, 1994, 2002, 2004, 2009]
#modoki_years = [1923, 1929, 1940, 1946, 1958, 1963, 1977]
modoki = []
for i in modoki_years:
    modoki.append(np.where(time==i)[0][0])

# El Nino years that were not modoki.
#elnino_years = [1963, 1965, 1972, 1982, 1987, 1994, 1997]
# All El nino years
elnino_years = [1911, 1913, 1914, 1918, 1925, 1930, 1941, 1951,
                1957, 1965, 1969, 1972, 1976, 1982, 1987, 1997, 2006]
#elnino_years = [1982, 1987, 1997, 2006]
#elnino_years = [1911, 1913, 1914, 1918, 1925, 1930, 1941, 1951,
#                1957, 1965, 1969, 1972, 1976]
elnino = []
for i in elnino_years:
    elnino.append(np.where(time==i)[0][0])


def composite_anomaly(data, index):
    """Composite data for index on first axis.
    """
    data_select = data[index,...].mean(axis=0)
    # Use 60:89 for 20cr and 50:79 for awap 1961-1990
    data_clim = data[base,...].mean(axis=0)
    return data_select - data_clim


# Composite modoki
hwf_modoki = composite_anomaly(hwf, modoki) 
hwd_modoki = composite_anomaly(hwd, modoki)
hwn_modoki = composite_anomaly(hwn, modoki)
hwm_modoki = composite_anomaly(hwm, modoki)
hwa_modoki = composite_anomaly(hwa, modoki)
hwt_modoki = composite_anomaly(hwt, modoki)

# Composite El Nino
hwf_elnino = composite_anomaly(hwf, elnino) 
hwd_elnino = composite_anomaly(hwd, elnino)
hwn_elnino = composite_anomaly(hwn, elnino)
hwm_elnino = composite_anomaly(hwm, elnino)
hwa_elnino = composite_anomaly(hwa, elnino)
hwt_elnino = composite_anomaly(hwt, elnino)

# Plot
def plotaus(data, signif, filename, units, title, rng):
    plt.figure()
    m = Basemap(projection='mill',
            llcrnrlon=110.,llcrnrlat=-45.,
            urcrnrlon=157.,urcrnrlat=-10.,
            resolution='l')
    lns, lts = np.meshgrid(lons, lats)
    x,y = m(lns,lts)
    cont = m.contour(x, y, data, linewidths=0.4, colors='k', levels=rng)
    for i,c in enumerate(cont.collections):
        if cont.levels[i]<0:
            c.set_dashes([(0, (2.0,2.0))])
    #plt.clabel(cont, fmt='%1.0f', fontsize=10)
    colors = m.pcolormesh(x, y, data, cmap='bwr', vmin=-rng[-1], vmax=rng[-1])
    m.contourf(x, y, signif, 1, colors='none', hatches=[None,'xx'])
    cbar = m.colorbar(colors, location='bottom', pad = 0.25)
    cbar.set_label(units)
    cbar.set_ticks(rng)
    plt.title(title)
    m.drawcoastlines()
    #m.drawmeridians(np.arange(110,157+10,10.),labels=[1,0,0,1],linewidth=0.2,fontsize=12)
    #m.drawparallels(np.arange(-40,0,10.),labels=[1,0,0,1],linewidth=0.2,fontsize=12)
    plt.savefig(filename, format='eps')
    plt.close()
    #plt.show()


#ts, pv = stats.ttest_ind(hwf[elnino],hwf[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwf_elnino, sig, 'hwf_elnino.svg', 'Days', 'HWF EP El Ni\~no', range(-16,17,4))
#
#ts, pv = stats.ttest_ind(hwn[elnino],hwn[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwn_elnino, sig, 'hwn_elnino.eps', 'Events', 'HWN El Nino', range(-10,11,2))
#
#ts, pv = stats.ttest_ind(hwd[elnino],hwd[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwd_elnino, sig, 'hwd_elnino.eps', 'Days', 'HWD El Nino', range(-15,16,5))
#
#ts, pv = stats.ttest_ind(hwa[elnino],hwa[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwa_elnino, sig, 'hwa_elnino.eps', 'degC^2', 'HWA El Nino', range(-20,21,5))
#
#ts, pv = stats.ttest_ind(hwm[elnino],hwm[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwm_elnino, sig, 'hwm_elnino.eps', 'degC^2', 'HWM El Nino', range(-20,21,5))
#
#ts, pv = stats.ttest_ind(hwf[modoki],hwf[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwf_modoki, sig, 'hwf_modoki.eps', 'Days', 'HWF Modoki', range(-15,16,3))
#
#ts, pv = stats.ttest_ind(hwn[modoki],hwn[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwn_modoki, sig, 'hwn_modoki.eps', 'Events', 'HWN Modoki', range(-10,11,2))
#
#ts, pv = stats.ttest_ind(hwd[modoki],hwd[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwd_modoki, sig, 'hwd_modoki.eps', 'Days', 'HWD Modoki', range(-15,16,5))
#
#ts, pv = stats.ttest_ind(hwa[modoki],hwa[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwa_modoki, sig, 'hwa_modoki.eps', 'degC^2', 'HWA Modoki', range(-20,21,5))
#
#ts, pv = stats.ttest_ind(hwm[modoki],hwm[50:79,...],axis=0,equal_var=False, nan_policy='omit')
#sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
#plotaus(hwm_modoki, sig, 'hwm_modoki.eps', 'degC^2', 'HWM Modoki', range(-20,21,5))


matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
#GRL Figures
pv = np.zeros(hwf.shape[1:])
for jj in xrange(hwf.shape[1]):
    for ii in xrange(hwf.shape[2]):
        _, pv[jj,ii] = stats.ks_2samp(hwf[elnino,jj,ii],hwf[base,jj,ii])
#ts, pv = stats.ttest_ind(hwf[elnino],hwf[50:79,...],axis=0,equal_var=False, nan_policy='omit')
sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
rng = range(-15,16,3)
f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,figsize=(6,10))
m = Basemap(ax=ax1, projection='mill',llcrnrlon=110.,llcrnrlat=-45.,urcrnrlon=157.,urcrnrlat=-10.,resolution='l',fix_aspect=True)
lns, lts = np.meshgrid(lons, lats)
x,y = m(lns,lts)
xx,yy = m(lns-.5,lts-.5)
cont = m.contour(x, y, hwf_elnino, linewidths=0.4, colors='k', levels=rng)
for i,c in enumerate(cont.collections):
    if cont.levels[i]<0:
        c.set_dashes([(0, (2.0,2.0))])
#plt.clabel(cont, fmt='%1.0f', fontsize=10)
colors = m.pcolormesh(xx, yy, hwf_elnino, cmap='bwr', vmin=-rng[-1], vmax=rng[-1])
m.contourf(x, y, sig, 1, colors='none', hatches=[None,'xxx'])
m.drawcoastlines()
m.drawmeridians(np.arange(110,157+10,10.),labels=[1,0,0,1],linewidth=0,fontsize=10)
m.drawparallels(np.arange(-40,0,10.),labels=[1,0,0,1],linewidth=0,fontsize=10)
cbar = m.colorbar(colors, location='bottom', pad=0.25)
cbar.set_label('Days')
cbar.set_ticks(rng)
for t in cbar.ax.get_xticklabels(): t.set_fontsize(10) 
ax1.set_title('a)', loc='left')
pv = np.zeros(hwf.shape[1:])
for jj in xrange(hwf.shape[1]):
    for ii in xrange(hwf.shape[2]):
        _, pv[jj,ii] = stats.ks_2samp(hwf[modoki,jj,ii],hwf[base,jj,ii])
sig = np.ma.array(pv<0.05, mask=hwf.mask[0])
m = Basemap(ax=ax2, projection='mill',llcrnrlon=110.,llcrnrlat=-45.,urcrnrlon=157.,urcrnrlat=-10.,resolution='l',fix_aspect=True)
lns, lts = np.meshgrid(lons, lats)
x,y = m(lns,lts)
m.drawmeridians(np.arange(110,157+10,10.),labels=[1,0,0,1],linewidth=0,fontsize=10)
m.drawparallels(np.arange(-40,0,10.),labels=[1,0,0,1],linewidth=0,fontsize=10)
cont = m.contour(x, y, hwf_modoki, linewidths=0.4, colors='k', levels=rng)
for i,c in enumerate(cont.collections):
    if cont.levels[i]<0:
        c.set_dashes([(0, (2.0,2.0))])
#plt.clabel(cont, fmt='%1.0f', fontsize=10)
colors = m.pcolormesh(xx, yy, hwf_modoki, cmap='bwr', vmin=-rng[-1], vmax=rng[-1])
m.contourf(x, y, sig, 1, colors='none', hatches=[None,'xxx'])
m.drawcoastlines()
cbar = m.colorbar(colors, location='bottom', pad=0.25)
cbar.set_label('Days')
cbar.set_ticks(rng)
for t in cbar.ax.get_xticklabels(): t.set_fontsize(10) 
ax2.set_title('b)', loc='left')
plt.savefig(filename='GRL_supfigure4.eps', format='eps')

# Calculate the area average diffference between modoki years and EP years.
neanc = nc.Dataset('/home/nfs/z5032520/analysis/model/mdknino_HWF_masks_awap.nc', 'r')
neaawap = neanc.variables['masks'][0,...].astype('bool')
neaawap = np.flipud(neaawap)
diff = hwf_elnino - hwf_modoki
neave = diff[neaawap].mean()
print 'Nino-Modoki Average over region:', neave
diff = hwf_elnino
neave = diff[neaawap].mean()
print 'Nino-Clim Average over region:', neave
diff = hwf_modoki
neave = diff[neaawap].mean()
print 'Modoki-Clim Average over region:', neave


# Plot all modoki years
#for i in modoki_years:
#    plotaus(hwf[np.where(time==i)[0][0]], sig, str(i)+'_heatwaves.svg', "Days", "HWF Modoki", range(-30,31,5))







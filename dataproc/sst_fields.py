import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc


# Load
pacninonc = nc.Dataset('pac_nino.nc')
pacnino = pacninonc.variables['sst'][:]
lats = pacninonc.variables['lat'][:]
lons = pacninonc.variables['lon'][:]
pacninanc = nc.Dataset('pac_nina.nc')
pacnina = pacninanc.variables['sst'][:]
indpiodnc = nc.Dataset('ind_piod.nc')
indpiod = indpiodnc.variables['sst'][:]
indniodnc = nc.Dataset('ind_niod.nc')
indniod = indniodnc.variables['sst'][:]

# Select a time series pont for each basin
i, j = 135, 72 # 253.125E, 0N
print 'lon: ',lons[i], 'lat: ', lats[j]
pacnino_se = pacnino[:,j,i]
pacnina_se = pacnina[:,j,i]
i, j = 41, 56 # 76.85E 20S
print 'lon: ',lons[i], 'lat: ', lats[j]
indpiod_se = indpiod[:,j,i]
indniod_se = indniod[:,j,i]

# set some global variables
units = '${}^{\circ}C$'
levels = np.arange(-2-0.1,2+0.2,0.2)

# plot figure
fig, axes = plt.subplots(nrows=3,ncols=2,figsize=(8.5,8))

# Pacific elnino
data = pacnino[11,...]
mp = Basemap(ax=axes[0,0],projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
cont = mp.contourf(xx,yy,data.squeeze(),
                   levels=levels,cmap='seismic',extend='both')
mp.drawcoastlines()
mp.drawparallels([-60,-30,0,30,60],labels=[True,False,False,True],linewidth=0.3)
mp.drawmeridians([30,180,330],labels=[True,False,False,True],linewidth=0.3)
plt.sca(axes[0,0])
plt.title('a)',loc='left')

# Pacific La Nina
data = pacnina[11,...]
mp = Basemap(ax=axes[1,0],projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
cont = mp.contourf(xx,yy,data.squeeze(),
                   levels=levels,cmap='seismic',extend='both')
mp.drawcoastlines()
plt.sca(axes[1,0])
plt.title('c)',loc='left')
cb = mp.colorbar(cont, location='bottom',pad=0.2)
cb.set_label(units)
cb.set_ticks(levels)
labels = ['-2.1','','','-1.5','','','-0.9','','','-0.3','','','0.3','','','0.9','','','1.5','','','2.1']
cb.set_ticklabels(labels)
mp.drawparallels([-60,-30,0,30,60],labels=[True,False,False,True],linewidth=0.3)
mp.drawmeridians([30,180,330],labels=[True,False,False,True],linewidth=0.3)

# Indian PIOD
data = indpiod[11,...]
mp = Basemap(ax=axes[0,1],projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
cont = mp.contourf(xx,yy,data.squeeze(),
                   levels=levels,cmap='seismic',extend='both')
mp.drawcoastlines()
mp.drawparallels([-60,-30,0,30,60],labels=[True,False,False,True],linewidth=0.3)
mp.drawmeridians([30,180,330],labels=[True,False,False,True],linewidth=0.3)
plt.sca(axes[0,1])
plt.title('b)',loc='left')

# Indian Niod
data = indniod[11,...]
mp = Basemap(ax=axes[1,1],projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
cont = mp.contourf(xx,yy,data.squeeze(),
                   levels=levels,cmap='seismic',extend='both')
mp.drawcoastlines()
mp.drawparallels([-60,-30,0,30,60],labels=[True,False,False,True],linewidth=0.3)
mp.drawmeridians([30,180,330],labels=[True,False,False,True],linewidth=0.3)
cb = mp.colorbar(cont, location='bottom',pad=0.2)
cb.set_label(units)
cb.set_ticks(levels)
labels = ['-2.1','','','-1.5','','','-0.9','','','-0.3','','','0.3','','','0.9','','','1.5','','','2.1']
cb.set_ticklabels(labels)
plt.sca(axes[1,1])
plt.title('d)',loc='left')

# Pacific Time series
plt.sca(axes[2,0])
plt.plot(pacnino_se,'k-',label='El Nino')
plt.scatter(np.arange(0,24), pacnino_se, marker='x',color='k')
plt.plot(pacnina_se,'k',label='La Nina', linestyle='dotted')
plt.scatter(np.arange(0,24), pacnina_se, marker='x',color='k')
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels.extend(labels)
plt.xticks(np.arange(0,24), labels, rotation=45,fontsize=9)
plt.axhline(y=0,color='k')
plt.xlim(0, 24)
plt.grid(True)
plt.xlabel('Month')
plt.ylabel('SST $^{\circ}$C')
plt.legend(loc='upper right',fontsize=9)
plt.sca(axes[2,0])
plt.title('e)', loc='left')

# Indian time series
plt.sca(axes[2,1])
plt.plot(indpiod_se,'k-',label='+IOD')
plt.scatter(np.arange(0,24), indpiod_se, marker='x',color='k')
plt.plot(indniod_se,'k',label='-IOD', linestyle='dotted')
plt.scatter(np.arange(0,24), indniod_se, marker='x',color='k')
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels.extend(labels)
plt.xticks(np.arange(0,24), labels, rotation=45,fontsize=9)
plt.axhline(y=0,color='k')
plt.xlim(0, 24)
plt.grid(True)
plt.xlabel('Month')
plt.ylabel('SST $^{\circ}$C')
plt.legend(loc='upper right',fontsize=9)
plt.sca(axes[2,1])
plt.title('f)', loc='left')

plt.tight_layout()

plt.savefig('sst_forcings.eps', format='eps')

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


# plot
def plot_map(data,lats,lons,units,rng):
    mp = Basemap(projection='robin', lon_0=180)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    mp.pcolormesh(xx,yy,data.squeeze(),
                  vmin=-rng,vmax=rng, cmap='seismic')
    cb = mp.colorbar(location='bottom')
    cb.set_label(units)
    mp.drawcoastlines()
    plt.show()


plot_map(pacnino[11,...],lats,lons,'${}^{\circ}C$',2)
plot_map(pacnina[11,...],lats,lons,'${}^{\circ}C$',2)
plot_map(indpiod[11,...],lats,lons,'${}^{\circ}C$',2)
plot_map(indniod[11,...],lats,lons,'${}^{\circ}C$',2)
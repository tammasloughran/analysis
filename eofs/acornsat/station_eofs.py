import __init__
import numpy as np
from netCDF4 import Dataset
from eofs.standard import Eof
from tools import rotate
from tools import pcaplot
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    # Load the data
    directory = '/srv/ccrc/data35/z5032520/ACORNSAT/'
    filename = 'EHF_heatwaves_ACORNSAT_1911_2014_summer.nc'
    ncin = Dataset(directory+filename, 'r')
    hwf = ncin.variables['HWF'][:]
    hwn = ncin.variables['HWN'][:]
    hwd = ncin.variables['HWD'][:]
    hwm = ncin.variables['HWM'][:]
    hwa = ncin.variables['HWA'][:]
    hwt = ncin.variables['HWT'][:]
    lat = ncin.variables['lat'][:]
    lon = ncin.variables['lon'][:]

    # Perform PCA
    retain = 4
    #metric_dict = {"HWN":hwn,"HWF":hwf,"HWD":hwd,"HWA":hwa,"HWM":hwm,"HWT":hwt}
    #for metric_name in metric_dict.keys():
    solver = Eof(hwn)
    pcs = solver.pcs(pcscaling=1, npcs=retain)
    explained_variance = solver.varianceFraction()
    errors = solver.northTest(vfscaled=True)
    eigens = solver.eigenvalues()
    eofs = solver.eofs(eofscaling=2, neofs=retain)
    eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
    eofs_correlation = solver.eofsAsCorrelation(neofs=retain)

    # Rotate. No need to reshape. 
    eofs2 = eofs
    pcs2 = pcs
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(eofs[0,...], eofs[1,...], eofs[2,...], marker='.')
    ax.plot([-1,1],[0,0],[0,0])
    ax.plot([0,0],[-1,1],[0,0])
    ax.plot([0,0],[0,0],[-1,1])
    ax.set_xlabel('EOF1')
    ax.set_ylabel('EOF2')
    ax.set_zlabel('EOF3')
    plt.title("Simple HWN EOFs")
    plt.show()
    eofs, R = rotate.varimax(eofs.T, normalize=False)
    pcs = np.dot(pcs, R)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(eofs.T[0,...], eofs.T[1,...], eofs.T[2,...], marker='.')
    ax.plot([-1,1],[0,0],[0,0])
    ax.plot([0,0],[-1,1],[0,0])
    ax.plot([0,0],[0,0],[-1,1])
    ax.set_xlabel('EOF1')
    ax.set_ylabel('EOF2')
    ax.set_zlabel('EOF3')
    plt.title("Rotated HWN EOFs")
    plt.show()

    # Make picures
    mapping = Basemap(projection='cyl',llcrnrlat=-44,
            urcrnrlat=-10,llcrnrlon=112, urcrnrlon=156)
    sign = np.sign(eofs2.T[:,0])
    mapping.scatter(lon,lat,s=abs(eofs2.T[:,0])*100,c=sign,cmap='bwr')
    mapping.drawcoastlines()
    plt.title("HWN simple EOF1")
    plt.savefig("HWN_simple_EOF1.eps",format='eps')
    fig = plt.figure()
    plt.plot(range(1911,2015), pcs2[:,0])
    plt.title("HWN simple PC1")
    plt.xlabel("Year")
    plt.ylabel("Normalized Score")
    plt.savefig("HWN_simple_PC1.eps",format='eps')

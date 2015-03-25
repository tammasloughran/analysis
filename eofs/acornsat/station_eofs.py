import __init__
import numpy as np
from netCDF4 import Dataset
from eofs.standard import Eof
from tools import rotate
from tools import pcaplot
import matplotlib.pyplot as plt

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
    eofs, R = rotate.varimax(eofs.T, normalize=False)
    pcs = np.dot(pcs, R)

    # Make picures
    plt.scatter(lon,lat,s=eofs[:,0]*70)
    plt.show()

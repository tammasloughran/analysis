from pyclimate.bpcca import BPCCA
from netCDF4 import Dataset
import matplotlib.pyplot as plt
tmaxfile = ('/media/Jupiter/reanalysis/AWAP/daily/tmax/tmax_detrended.nc')
hwfile = ('/media/Jupiter/reanalysis/AWAP/yearly/ehfhw/CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
hwnc = Dataset(hwfile,'r')
tmaxnc = Dataset(tmaxfile,'r')
tmax = tmaxnc.variables['tmax'][-104:]
hw = hwnc.variables['HWF_EHF'][:]
CCA = BPCCA(tmax, hw, (4,4))
L = CCA.leftPatterns()
R = CCA.rightPatterns()
r = CCA.correlation()
a = CCA.rightExpCoeffs()
b = CCA.leftExpCoeffs()
plt.contourf(R[:,:,0])
plt.show()
plt.plot(a)
plt.show()

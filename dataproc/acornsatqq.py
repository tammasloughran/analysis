from netCDF4 import Dataset
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

# Load data.
ncin = Dataset('/srv/ccrc/data35/z5032520/ACORNSAT/daily/ACORNSAT.nc','r')
tmax = ncin.variables['tmax'][:]
tmin = ncin.variables['tmin'][:]
stationid = ncin.variables['stationid'][:]

# Calculate percentage of existing data over the entire period.
tmax_percent_exist = (1-(np.isnan(tmax).sum(axis=0)/float(len(tmax))))*100.
tmin_percent_exist = (1-(np.isnan(tmin).sum(axis=0)/float(len(tmin))))*100.

# Test if stations have >= 80% of data over the analysis period.
tmax_pass_test = (tmax_percent_exist >= 80.)
tmin_pass_test = (tmin_percent_exist >= 80.)

# Do tmax and tmin agree on bad stations?
agreement_test = (tmax_pass_test==tmin_pass_test).all()
if agreement_test:
    print "The following stations have < 80% of data vailable:"
    print stationid[tmax_pass_test==False]
    np.save('qqtest.npy',tmax_pass_test)
    print "qqtest.npy file saved"
else:
    print "Tmax and Tmin have different stations that have < 80% of data."
    print "For Tmax they are:"
    print stationid[tmax_pass_test==False]
    print "For Tmin they are:"
    print stationid[tmin_pass_test==False]

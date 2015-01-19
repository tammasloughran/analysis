'''Calculate the EOFs of EHF heatwave metrics.

EHF heat wave metrics have been calculated from AWAP Tmax and Tmin data using a
90th percentile with respect to the entire analysis period. Metrics are annual
heat wave characteristics. They are HWF (frequency), HWD (duration), HWA 
(amplitude), HWM (magnitude), HWN (number) and HWT (timing).
'''

def load_heat_waves(filename):
    '''Load heat wave metrics from a netcdf file.

    Arguments 
    filename -- name of file containing heat wave metrics.
    maskname -- name of file containing AWAP land-sea mask.

    Returns
    hwf -- frequency
    hwn -- number
    hwd -- duration
    hwa -- amplitude
    hwm -- magnitude
    hwt -- timing
    '''
    from netCDF4 import Dataset
    ncin = Dataset(filename, 'r')
    hwf = ncin.variables['HWF_EHF'][:]
    hwn = ncin.variables['HWN_EHF'][:]
    hwd = ncin.variables['HWD_EHF'][:]
    hwa = ncin.variables['HWA_EHF'][:]
    hwm = ncin.variables['HWM_EHF'][:]
    hwt = ncin.variables['HWT_EHF'][:]
    lat = ncin.variables['lat'][:]
    lon = ncin.variables['lon'][:]
    times = ncin.variables['Times'][:]
    return hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times

if __name__ == "__main__":
    from numpy import arange, cos, sqrt, deg2rad, newaxis, dot
    from eofs.standard import Eof
    import rotate
    import pcaplot
    start_year = 1911
    end_year = 2014
    # Load the heat wave metrics.
    fname = ("/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/CCRC_NARCliM_1911-"
             "2014_EHFheatwaves_summer_AWAP0.5deg.nc")
    hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times = load_heat_waves(fname)
    # Calculate weightings.
    coslat = cos(deg2rad(lat)).clip(0.,1.)
    wgts = sqrt(coslat)[..., newaxis]
    # Set up solver
    retain = 4
    solver = Eof(hwf, weights=wgts)
    pcs = solver.pcs(pcscaling=1, npcs=retain)
    explained_variance = solver.varianceFraction()
    errors = solver.northTest(vfscaled=True)
    eigens = solver.eigenvalues()
    eofs = solver.eofs(eofscaling=2, neofs=retain)
    eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
    eofs_correlation = solver.eofsAsCorrelation(neofs=retain)
    # Apply rotation to PCs and EOFs.
    eofs2 = eofs
    pcs, eofs = rotate.do_rotation(pcs, eofs, space='state')
    # Plotting.
    years = arange(start_year,end_year+1)
    pcaplot.plot_eigenvalues(explained_variance, errors)
    pcaplot.plot_eofs(eofs, lon, lat, 'Rotated_EOFs')
    pcaplot.plot_eofs(eofs2, lon, lat, 'EOFs')
    pcaplot.plot_eofs(eofs_covariance, lon, lat, 'EOFs_Covariance')
    pcaplot.plot_eofs(eofs_correlation, lon, lat, 'EOFs_Correlation')
    pcaplot.plot_pcs(pcs, years)
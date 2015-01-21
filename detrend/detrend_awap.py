'''Detrend AWAP tmax and tmin.
'''

def loadnc(infile,varname):
    '''Load the netCDF file and return varname data.
    '''
    from netCDF4 import Dataset
    ncfile = Dataset(infile, 'r')
    data = ncfile.variables[varname][:]
    return ncfile, data

def copync(outfile, copyfile, exclude=''):
    '''Copy a netCDF file structure.

    Arguments
    outfile -- filename string of new netCDF file.
    copyfile -- netCDF4 Dataset object to be copied.
    exclude -- a list of variables that should be excluded from copy.

    Returns
    outnc -- a netCDF4 Dataset object copy.
    '''
    from netCDF4 import Dataset
    outnc = Dataset(outfile, 'w')
    # Copy attributes
    for attribute in copyfile.ncattrs():
        outnc.setncattr(attribute, copyfile.getncattr(attribute))
    # Copy dimensions
    for dimname, dim in copyfile.dimensions.iteritems():
        outnc.createDimension(dimname, len(dim))
    # Copy variables
    for varname, var in copyfile.variables.iteritems():
        copyvar = outnc.createVariable(varname, var.datatype, var.dimensions)
        # Copy variable attributes
        for attribute in copyfile.variables[varname].ncattrs():
            copyvar.setncattr(attribute, \
                    copyfile.variables[varname].getncattr(attribute))
        # Exclude variables in exclude
        if not any(varname in string for string in exclude):
            copyvar[:] = var[:]
    return outnc

from scipy.signal import detrend
# Load data
tmaxfile = ('/srv/ccrc/data35/z5032520/AWAP/daily/tmax/'
           'AWAP_TX_1911-2014_0.5deg.nc')
tminfile = ('/srv/ccrc/data35/z5032520/AWAP/daily/tmin/'
           'AWAP_TN_1911-2014_0.5deg.nc')
tmaxnc, tmax = loadnc(tmaxfile, 'tmax')
tmax_detrended = detrend(tmax)
outfilename = 'tmax_detrended.nc'
tmaxoutnc = copync(outfilename, tmaxnc, exclude='tmax')
tmaxoutvar = tmaxoutnc.variables['tmax']
tmaxoutvar[:] = tmax_detrended[:]
tminnc, tmin = loadnc(tminfile, 'tmin')
tmin_detrended = detrend(tmin)
outfilename = 'tmin_detrended.nc'
tminoutnc = copync(outfilename, tminnc, exclude='tmin')
tminoutvar = tminoutnc.variables['tmin']
tminoutvar[:] = tmin_detrended[:]


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

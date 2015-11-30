def bh_fdr(p,a):
    """Determine which of unsorted p values reject the null hypothesis using 
    the Benjamini & Hoch false discovery rate.
    
    Arguments:
    p - an array of p values.
    a - the global significance threshold.

    Returns:
    sig_unsort - unsorted indicator of significance.
    p_fdr - the maximum p value for the false discovery rate.
    """
    import numpy as np
    # Flatten the array and store the original shape
    oshape = p.shape
    if len(p.shape)>1:
        p = p.flatten()
    # If there are nans then exlcule these
    hasnan = False
    if np.isnan(p).any():
        hasnan = True
        nanshape = p.shape
        nan_data = np.isnan(p)
        p = p[np.logical_not(nan_data)]
    # Get the number of tests
    n = len(p)
    # Sort the p values
    ps = np.sort(p)
    j = np.arange(1,n+1)
    # Calculate the critical FDR rates for each p value
    crit = (j/float(n))*a
    # determine significance for the sorted p values
    sig = ps<=crit
    # find the maximum p value for the significant FDR
    if np.logical_not(sig).all():
        p_fdr = None
    else:
        p_fdr = np.max(ps[sig])
    # Unsort the significance values back into the original order
    sig_unsort = np.ones(len(sig), dtype=bool)
    for ii, pv in enumerate(ps):
        i = np.where(p==pv)[0][0]
        sig_unsort[i] = sig[ii]
    # Inset nans back into the data.
    if hasnan:
        dummy = np.ones(nanshape)*np.nan
        dummy[np.logical_not(nan_data)] = sig_unsort
        sig_unsort = dummy
    # Reshape significance back into input dimensions.
    if oshape!=0: sig_unsort = sig_unsort.reshape(oshape)
    return sig_unsort, p_fdr

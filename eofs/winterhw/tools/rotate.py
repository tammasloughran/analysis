"""
This module provides functions for performing varmiax rotations on
emperical orthogaonal functions. Tammas F. Loughran
"""

def do_rotation(pcs, eofs):
    '''Prepare data and perform varimax rotation on EOF loadings.

    The returned rotated PC scores are only meaningful if the simple EOFs 
    provided are loadings that have been scaled by multiplying each one by 
    the square root of it's respective eigenvalue.
    
    The EOFs may have missing values. These are removed prior to rotation and
    returned after.

    Arguments
    pcs -- matrix of PCs allong the second axis.
    eofs -- matrix of EOFs allong the first axis.

    Returns
    pcs -- rotated PCs.
    eofs -- rotated EOFs.
    '''
    from numpy import dot, ma, where, ones, isnan, NaN
    import numpy as np
    # Reshape
    if len(eofs.shape)>2:
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx * ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
    else:
        eofs2d = eofs
    # Remove missing
    nonMissingIndex = where(isnan(eofs2d.data[0]) == False)[0]
    eofs2d = eofs2d.data[:, nonMissingIndex]
    # Rotate
    rot_eofs_nomiss, R = varimax(eofs2d.T, normalize=False)
    # Restore missing
    rotated_eofs = ones([nmaps, ngridpoints]) * NaN
    rotated_eofs = rotated_eofs.astype(eofs2d.dtype)
    rotated_eofs[:, nonMissingIndex] = rot_eofs_nomiss.T
    # Reshape
    if len(eofs.shape)>2:
        rotated_eofs = rotated_eofs.reshape([nmaps, ny, nx])
    # Restore
    eofs = ma.masked_array(rotated_eofs, eofs.mask)
    # Calculate rotated PCs
    pcs = dot(pcs-np.mean(pcs,axis=0), R)
    return pcs, eofs, R

def expvar_from_eofs(eofs,totalvar):
    """Calculate the explained variance of a set of eofs compared the
    total explained variance.
    """
    import numpy as np
    # Reshape
    if len(eofs.shape)>2: 
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx * ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
    else:
        eofs2d = eofs
    # Remove missing
    nonMissingIndex = np.where(np.isnan(eofs2d.data[0]) == False)[0]
    eofsNoMissing = eofs2d.data[:, nonMissingIndex]
    # Calculate explained variance
    expvar = np.sum(np.power(eofsNoMissing,2),axis=1)/totalvar
    return expvar

def project_eofs(eofs, field, totalvar):
    """Project eof loadings with missing values onto a field.
    Only use this if the rotation was performed on raw EOF loadings.
    (i.e. the eigenvectors that have no scaling.)
    """
    import numpy as np
    # Get dimensions
    samples = field.shape[0]
    nmaps, ny, nx = eofs.shape
    ngridpoints = nx * ny
    # Reshape to 2 dimensions and remove missing values.
    eofs2d = eofs.reshape([nmaps, ngridpoints])
    nonMissingIndex = np.where(np.isnan(eofs2d.data[0]) == False)[0]
    eofsNoMissing = eofs2d.data[:, nonMissingIndex]
    field2d = field.reshape([samples, ngridpoints])
    fnonMissingIndex = np.where(np.isnan(field2d.data[0]) == False)[0]
    fNoMissing = field2d.data[:, nonMissingIndex]
    # Calculate PCs
    pcs = np.dot(fNoMissing,eofsNoMissing.T)
    return pcs

def varimax(X, normalize=True, gamma = 1.0, it = 200, tol = 1e-7):
    """Performs a varimax rotation on the input matrix.
    
    The method used uses a singular value decomposition (SVD) of 
    the varimax criterion. The rotation matrix is found itteratively.

    Arguments
    X     - A pxk input matrix that will be rotated.
    normalize - If True normalize matrix for rotation & denormalise after.
    gamma - Specifies the rotation type. 1 for varimax and 0 for quartimax.
    it    - Maximum number of itterations.
    tol   - Tollerance of convergence.

    Local variables.
    p,k   - Dimensions of the input matrix
    R     - Rotation matrix that we want to find. The initial 
            estimation for R is R = I_k (i.e. an identity matrix that
            rotates by nothing.)
    sums  - A measure to test if the solution has converged. 
    Lambda - Itteratively rotated input matrix.
    D     - Diagonal matrix containing the sum of the squared weightings.
    C     - Cubed weightings.
    U,S,VT - Output of the SVD. S is the scalings & U*VT give the rotation.
    
    This is a slight modification if the R implementation if the 
    varimax function. The goal of this python implementation is to produce
    well documented code to aid clarity and understanding of varimax
    rotation. The matrices C and D are consistent with the explanation of
    the varimax method given in chapter 6 of "Factor Analysis as a
    Statistacal Method" D. N. Lawley & A. E. Maxwell, 1963, London
    Butterworths. The varimax rotation criterion is

    A = L'*C - (1/p)*L'*L*D

    Example
    B, R = rotate.varimax(N)
    """
    from numpy import eye, dot, sum, diag, sqrt, newaxis
    from numpy.linalg import svd

    # Normalize the matrix.
    if normalize==True:
        h = sqrt(sum(X**2,1))
        X = X / h[:,newaxis]

    # Rotate the matrix.
    p,k = X.shape
    R = eye(k)
    sums=0
    for i in xrange(it):
        sums_old = sums
        Lambda = dot(X, R)
        D = diag(sum(Lambda**2,0))
        C = Lambda**3
        U,S,VT = svd(dot(X.T, C - (gamma/p) * dot(Lambda, D)))
        R = dot(U,VT)
        sums = sum(S)
        if sums_old!=0 and sums/sums_old < 1 + tol: break
    X = dot(X,R)

    # Denormalize.
    if normalize==True:
        X = X * h[:,newaxis]

    return X, R

"""
This module provides functions for performing varmiax rotations on
emperical orthogaonal functions or principle components.
"""

def do_rotation(pcs, eofs, space='State'):
    '''Prepare data and perform varimax rotation.
    
    do_rotation reshapes the EOFs for rotation, then applies the 
    rotation and reshapes back into the original form.
    Rotation can be done on the PCs in sample space or on the EOFs in
    state space. In state space the masked values in the EOFs are removed
    before rotation and returned afterwards.

    Arguments
    pcs -- matrix of PCs allong the second axis.
    eofs -- matrix of EOFs allong the first axis.
    space -- type of varimax rotaion. Either 'state' or 'sample'.

    Returns
    pcs -- rotated PCs.
    eofs -- rotated EOFs.
    '''
    from numpy import dot, ma, where, ones, isnan, NaN
    if space == 'sample':
        pcs, R = varimax(pcs)
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx*ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
        rot_eofs = dot(R, eofs2d)
        eofs = rot_eofs.reshape([nmaps, ny, nx])
    elif space == 'state':
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx * ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
        nonMissingIndex = where(isnan(eofs2d.data[0]) == False)[0]
        dataNoMissing = eofs2d.data[:, nonMissingIndex]
        rot_eofs_nomiss, R = varimax(dataNoMissing.T, normalize=False)
        rotated_eofs = ones([nmaps, ngridpoints]) * NaN
        rotated_eofs = rotated_eofs.astype(eofs2d.dtype)
        rotated_eofs[:, nonMissingIndex] = rot_eofs_nomiss.T
        rotated_eofs = rotated_eofs.reshape([nmaps, ny, nx])
        eofs = ma.masked_array(rotated_eofs, eofs.mask)
        pcs = dot(pcs, R)
    return pcs, eofs

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
    B = rotate.varimax(N)

    Tammas F. Loughran
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

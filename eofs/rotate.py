def varimax(X, gamma = 1.0, q = 20, tol = 1e-6):
    """
    varimax() performs an varimax rotation on the input matrix 
    and returns the rotated matrix. The method used uses a 
    singular value decomposition (SVD) of the varimax criterion. The
    rotation matrix is found itteratively.

    Inputs
    X     - A pxk input matrix that will be rotated.
    gamma - Specifies the rotation type. 1 for varimax and 0 for quartimax.
    q     - Maximum number of itterations.
    tol   - Tollerance of convergence.

    Local variables.
    p,k   - Dimensions of the input matrix
    R     - Rotation matrix that we want to find. The initial 
            estimation for R is R = I_k (i.e. an identity matrix that
            rotates by nothing.)
    sums  - A measure to test if the solution has converged. 
    s     - Contains the singular values of the varimax criterion SVD.
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
    Butterworths.

    Tammas F. Loughran
    """
    from numpy import eye, asarray, dot, sum, diag
    from numpy.linalg import svd
    p,k = X.shape
    R = eye(k)
    sums=0
    for i in xrange(q):
        sums_old = sums
        Lambda = dot(X, R)
        D = diag(sum(Lambda**2,0))
        C = Lambda**3
        U,S,VT = svd(dot(X.T, C - (gamma/p) * dot(Lambda, D)))
        R = dot(U,VT)
        sums = sum(S)
        if sums_old!=0 and sums/sums_old < 1 + tol: break
    return dot(X, R)

def kaiser_normalise(M):
    """kaiser_normalise() scales the weightings of the input matrix according to the kaiser normalisation.
    Kaiser normaisation divides each element of the matrix by h_i, where h_i^2 is the communality if the jth test.
    The communality is simply the sum of the square of the weightings for a given pc.

    DO NOT USE THIS YET! ITS COMPLETELY WRONG.

    """
    from numpy import array
    import math
    mshape = array(M.shape)
    for j in range(mshape[1])
        pc = M[:,j]
        h_i2 = sum(pc**2)
        M[:,j] = M[:,j] / math.sqrt(h_i2)
    return M 

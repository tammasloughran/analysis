"""
This module provides functions for performing varmiax rotations on
principle components or factor analasis.
"""

def varimax(X, normalize=True, gamma = 1.0, it = 200, tol = 1e-7):
    """
    varimax(X) performs a varimax rotation on the input matrix
    and returns the rotated matrix. The method used uses a
    singular value decomposition (SVD) of the varimax criterion. The
    rotation matrix is found itteratively.

    Inputs
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
    from numpy import eye, dot, sum, diag
    from numpy.linalg import svd
    # Normalize the matrix.
    if normalize==True:
        X, h_i = kaiser_normalize(X)
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
        X = kaiser_denormalize(X, h_i)
    return X

def kaiser_normalize(X):
    """
    kaiser_normalize(X) scales the weightings of the input matrix 
    according to the kaiser normalization. 

    Input
    X - The matrix to be normalized.

    Output
    X - The normalized matrix.
    h - A vector containing the communalities of the unscaled X.

    Kaiser normaization divides each row of the matrix by h_i, 
    where h_i^2 is the communality. The communality is simply 
    the sum of the square of the weightings for each row.

    Example
    N,h1 = rotate.kaiser_normalize(A)
    """
    from numpy import sum, sqrt, newaxis
    h = sqrt(sum(X**2,1))
    X = X / h[:,newaxis]
    return X, h

def kaiser_denormalize(X, h):
    """
    kaiser_denormalize(X, h) denormalizes a matrix that has been rotated.
    It arguments are a matix to denormalise and a vector of 
    communalities. The length of the vector containing the communalities 
    must match with the number of rows of the input matrix.

    Input
    X - Matrix to denormalize
    h - communality vector from kaiser normalize.

    Output
    X - denormalized matrix.

    Example
    Ar = rotate.kaiser_denormalize(B,h1)
    """
    from numpy import newaxis
    class CommunalitiesMismatch(Exception):
        """
        This exception returns the length of h and the length of
        second axis of X.
        """
        def __init__(self, lenh, xshape):
            self.shapes = [lenh, xshape]
        def __str__(self):
            return repr(self.shapes)
    pass
 
    # Check if the lengths agree.
    if len(h) != X.shape[0]:
        raise CommunalitiesMismatch(len(h), X.shape[0])
    
    # Denormalize
    X = X * h[:,newaxis]
    return X

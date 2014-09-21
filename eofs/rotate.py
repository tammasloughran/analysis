def varimax(Phi, gamma = 1.0, q = 20, tol = 1e-6):
    from numpy import eye, asarray, dot, sum 
    from numpy.linalg import svd
    p,k = Phi.shape
    R = eye(k)
    d=0
    for i in xrange(q):
        d_old = d
        Lambda = dot(Phi, R)
        u,s,vh = svd(dot(Phi.T,asarray(Lambda)**3 - (gamma/p) * dot(Lambda, diag(diag(dot(Lambda.T,Lambda))))))
        R = dot(u,vh)
        d = sum(s)
        if d_old!=0 and d/d_old < 1 + tol: break
    return dot(Phi, R)

def kaiser_norm(M):
    """normalise() scales the weightings of the input matrix according to the kaiser normalisation.
    Kaiser normaisation divides each element of the matrix by h_i, where h_i^2 is the communality if the jth test.
    The communality is simply the sum of the square of the weightings for a given pc."""
    from numpy import array
    import math
    mshape = array(M.shape)
    for j in range(mshape[1])
        pc = M[:,j]
        h_i2 = sum(pc**2)
        M[:,j] = M[:,j] / math.sqrt(h_i2)
    return M 

#!/usr/bin/python3.11

from random import random, seed
from numpy import array, transpose, allclose, conj, dot, sqrt
from copy import deepcopy

def SysMake(size, min, max):
    #Make coefficient matrix A
    mat = [[None for i in range(0, size)] for j in range(0, size)]

    for i in range(0, size):
        mat[i][i] = min + random()*(max-min)
        for j in range(i + 1, size):
            re = min + random()*(max-min)
            im = min + random()*(max-min)
            mat[i][j] = re + 1j*im
            mat[j][i] = re - 1j*im

    ##and check if it is really hermitian
    if(not allclose(array(mat), conj(transpose(mat)))):
        print("I suck")
        exit(1)

    #Make resultant vector b
    res = [None for i in range(0, size)]
    for i in range(0, size):
        re = min + random() * (max - min)
        im = min + random() * (max - min)
        res[i] = re + 1j*im

    return mat, res

def CAV(z): ##absolute value of a complex number
    return sqrt(z.real*z.real + z.imag*z.imag)

def IndTrns(n, i, j):
    return n*i + j



from numpy import array
#Based on https://github.com/fredrik-johansson/mpmath/blob/master/mpmath/matrices/eigen_symmetric.py
def QR(n, A, giveVecs=False, smallTol=1e-8, itrLim=30):

    #ints
    i, j, k, l, m = None, None, None, None, None
    iterCount = None

    #real floats
    sc, scinv, f1, h, g, hhparam = None, None, None, None, None, None
    rtmp, r, s, g, c, p, b, f = None, None, None, None, None, None, None, None
    #complex floats
    ctmp, f0 = None, None

    #real arrays
    d = [0. for i in range(0, n)]
    e = [0. for i in range(0, n)]
    #complex arrays
    t = [0. + 0j for i in range(0, n)]


    ##Transform 2D array into 1D
    a_ = [None for i in range(0, n*n)]
    for i in range(0, n):
        a_[IndTrns(n, i, i)] = A[i][i]
        for j in range(i + 1, n):
            #a_[IndTrns(n, i, j)] = conj(A[i][j])
            #a_[IndTrns(n, j, i)] = A[i][j]
            a_[IndTrns(n, i, j)] = A[i][j]
            a_[IndTrns(n, j, i)] = conj(A[i][j])
    A = array(a_)



    ##################################
    ### TRIDIAG POTION OF THE ALGO ###
    ##################################
    t[n-1] = 1. + 0.j
    for i in range(n-1, 0, -1):
        l = i - 1

        ##vector scaling
        sc = 0.
        for j in range(0, i):
            ij = IndTrns(n, i, j)
            sc += abs(A[ij].real) + abs(A[ij].imag)

        ###early term if scale is small - no need for useless, expensive calculations
        if(abs(sc) < smallTol):
            e[i] = 0.
            d[i] = 0.
            t[l] = 1. + 0j
            continue

        scinv = 1./sc


        ##stopping criteria for last iteration
        if(i == 1):
            f0 = A[IndTrns(n, i, l)]
            f1 = abs(f0)
            e[i] = f1
            d[i] = 0.
            if(f1 > smallTol):
                t[l] = t[i]*f0/f1
            else:
                t[l] = t[i]
            continue


        ##householder transformation setup
        h = 0.
        for j in range(0, i):
            ij = IndTrns(n, i, j)
            A[ij] *= scinv
            h += A[ij].real*A[ij].real + A[ij].imag*A[ij].imag

        f0 = A[IndTrns(n, i, l)]
        f1 = abs(f0) ##f1 = CAV(f0) ???
        g = sqrt(h)
        h += g*f1
        e[i] = sc*g
        if(f1 > smallTol):
            f0 /= f1
            ctmp = -t[i]*f0
            g *= f0
        else:
            ctmp = -t[i]
        A[IndTrns(n, i, l)] += g
        f0 = 0.


        ##householder transformation application
        for j in range(0, i):
            if(giveVecs):
                A[IndTrns(n, j, i)] = A[IndTrns(n, i, j)] / h

            ###A dot U
            g = 0 + 0j
            for k in range(0, j+1):
                g += conj(A[IndTrns(n, j, k)]) * A[IndTrns(n, i, k)]
            for k in range(j+1, i):
                g += A[IndTrns(n, k, j)] * A[IndTrns(n, i, k)]

            ###P
            t[j] = g/h
            f0 += conj(t[j])*A[IndTrns(n, i, j)]

        hhparam = 0.5*f0/h
        ###reduce A, get Q
        for j in range(0, i):
            f0 = A[IndTrns(n, i, j)]
            g = t[j] - hhparam*f0
            t[j] = g

            for k in range(0, j + 1):
                A[IndTrns(n, j, k)] -= conj(f0)*t[k] + conj(g)*A[IndTrns(n, i, k)]

        t[l] = ctmp
        d[i] = h
    #print(d)

    #################################
    ### EIGVAL POTION OF THE ALGO ###
    #################################
    def safehypot(a, b): ##aux function for (a^2 + b^2)^(1/2) w/o under/overflow
        #return sqrt(a*a + b*b)
        absa, absb = abs(a), abs(b)
        if(absa > absb):
            return absa*(1. + (absb/absa)**(2))**(1/2)
        if(absb == 0.0):
            return 0.0
        return absb*(1. + (absa/absb)**(2))**(1/2)

    #setup arrays: shift {e} to the left, make d explicitly real-valued
    for i in range(1, n):
        e[i-1] = e[i]
    e[n-1] = 0.

    d[0] = 0.
    for i in range(0, n):
        ii = IndTrns(n, i, i)
        rtmp = d[i]
        d[i] = A[ii].real
        A[ii] = rtmp
    #note that at this point, all values are strictly real

    B = [[0. + 0j for j in range(0, n)] for i in range(0, n)]
    for i in range(0, n):
        B[i][i] = 1. + 0j
    #from copy import deepcopy
    #B = deepcopy(A)
    B = array(B)

    for l in range(0, n):
        itrCount = 0


        while(1):
            ##grab a small off-diagonal element
            m = l
            while(1):
                if(m + 1 == n):
                    break
                if(  abs(e[m]) < smallTol*(abs(d[m]) + abs(d[m+1]))  ):
                    break
                m += 1
            if(m == l):
                break

            if(itrCount >= itrLim):
                print("I suck")
                exit(2)

            itrCount += 1

            ##shift
            g = 0.5*(d[l+1] - d[l])/e[l]
            r = safehypot(g, 1.0)
            s = g - r if g < 0.0 else g + r
            g = d[m] - d[l] + e[l] / s

            ##plane-into-givens rotations to get back to tridiagonal form
            s, c, p = 1., 1., 0.
            for i in range(m-1, l-1, -1):
                f = s*e[i]
                b = c*e[i]

                ###slightly better than the standard QL algo
                if(abs(f) > abs(g)):
                    c = g/f
                    r = safehypot(c, 1.)
                    e[i+1] = f*r
                    s = 1/r
                    c *= s
                else:
                    s = f/g
                    r = safehypot(s, 1.)
                    e[i+1] = g*r
                    c = 1/r
                    s *= c

                g = d[i+1] - p
                r = (d[i] - g)*s + 2.*c*b
                p = s*r
                d[i+1] = g + p
                g = c*r - b

                ##Get eigenvectors if necessary
                if(giveVecs):
                    for w in range(0, n):
                        f = B[w][i+1]
                        B[w][i+1] = s*B[w][i] + c*f
                        B[w][i]   = c*B[w][i] - s*f

            #if(abs(r) <= smallTol and i >= l):
            #    continue
            d[l] -= p
            e[l] = g
            e[m] = 0.

    #sorting: NOTE: BUBBLE SORT - DON'T USE THIS GARBAGE
#    for ii in range(1, n):
#        i = ii - 1
#        k = i
#        p = d[i]
#        for j in range(ii, n):
#            if(d[j] >= p):
#                continue
#            k = j
#            p = d[k]
#        if(k == i):
#            continue
#        d[k] = d[i]
#        d[i] = p
#
#        if(giveVecs):
#            for w in range(0, n):
#                p = B[w][i]
#                B[w][i] = B[w][k]
#                B[w][k] = p

    if(not giveVecs):
        return d ##d holds the (unsorted) eigenvalues
    
    for i in range(0, n):
        for k in range(0, n):
            B[k][i] *= t[k]
    for i in range(0, n):
        if(A[IndTrns(n, i, i)] != 0):
            for j in range(0, n):
                ctmp = 0 + 0j
                for k in range(0, i):
                    ctmp += conj(A[IndTrns(n, k, i)]) * B[k][j]
                    #print(conj(A[IndTrns(n, k, i)]), B[k][j], conj(A[IndTrns(n, k, i)]) * B[k][j], ctmp)
                    print(ctmp)
                for k in range(0, i):
                    B[k][j] -= ctmp*A[IndTrns(n, i, k)]

    exit(6)
    return array(d), conj(array(B)) ##just the conj, NOT transpose!

"""
**************************************************************************************************************
*****  TESTING ***********************************************************************************************
**************************************************************************************************************
"""

from numpy import round, dot, array, allclose
import time

SIZE = 3
MIN, MAX = -10, +10
seed(9642069)
A, b = SysMake(SIZE, MIN, MAX)
#A = [[1, 2 + 5j], [2 - 5j, 3]]
A, b = array(A), array(b)
print(A)

from numpy.linalg import eigh
e, v = eigh(A)
#print(allclose(dot(A, v) - e*v, 0.0 + 0.0j, 1e-6, 1e-6))
#print(dot(A, vecs) - e*array(vecs))
#print(round(v, 3))
#print(round(e, 1), "\n\n")

start = time.time()
eigs, vecs = QR(SIZE, A, True)
end = time.time()
print(allclose(dot(A, vecs) - eigs*array(vecs), 0.0 + 0.0j, 1e-6, 1e-6))
#print(dot(A, vecs) - e*array(vecs))
print(round(vecs, 3))

print(allclose(sorted(e), sorted(eigs), 1e-6, 1e-6))
#print(allclose(e, eigs, 1e-6, 1e-6))

print("TIME:", end-start)
#print(round(vecs, 2))


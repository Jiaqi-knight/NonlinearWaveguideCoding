# -*- coding: utf-8 -*-

'''
=====================================================================================

Copyright (c) 2017 - 2018 Université de Lorraine & Luleå tekniska universitet
Author: Luca Di Stasio <luca.distasio@gmail.com>
                       <luca.distasio@ingpec.eu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=====================================================================================

DESCRIPTION

Tested with Python 2.7 Anaconda 2.4.1 (64-bit) distribution in Windows 7.

'''

import sys

def zeros(m,n):
    res = []
    if m>1 and n>1:
        for i in range(0,m):
            row = []
            for j in range(0,n):
                row.append(0.0)
            res.append(row)
    elif m>1:
        for i in range(0,m):
            res.append(0.0)
    else:
        for i in range(0,n):
            res.append(0.0)
    return res

def zeroTensor(indexList):
    res = []
    if len(indexList)>1:
        element = zeroTensor(indexList[1:])
        for i in range(indexList[0]):
            res.append(element)
    else:
        for i in range(indexList[0]):
            res.append(0.0)
    return res

def ones(m,n):
    res = []
    if m>1 and n>1:
        for i in range(0,m):
            row = []
            for j in range(0,n):
                row.append(1.0)
            res.append(row)
    elif m>1:
        for i in range(0,m):
            res.append(1.0)
    else:
        for i in range(0,n):
            res.append(1.0)
    return res

def dividedDiffNewton(x,y):
    #  Input: N x 1 vector x of interpolation nodes
    #         N x 1 vector y of values at interpolation nodes
    #  Output: N x N lower-triangular matrix d of divided differences
    my = len(y)
    d = zeros(my,my);
    for i in range(0,my):
        d[i][0] = y[i]
    for j in range(1,my):
        for i in range(j,my):
            d[i][j] = (d[i][j-1]-d[i-1][j-1])/(x[i]-x[i-j+1])
    return d

def interpolantLagrange(x,xval,index):
    #  Input: N x 1 vector x of interpolation nodes
    #         M x 1 vector xval of nodes for function evaluation
    #         index "index" of current node
    #  Output: M x 1 vector l of interpolant evaluations
    N = len(x);
    M = len(xval);
    l = ones(M,1);
    for i in range(0,N):
        if i!=index:
            for j in range(0,M):
                l[j] = l[j]*(xval[j]-x[i])/(x[index]-x[i])
    return l

def interpLagrange1D(x,y,z):
    # Input: N x 1 vector x of interpolation nodes
    #        N x 1 vector y of values at interpolation nodes
    #        M x 1 vector z of nodes for function evaluation
    #  Output: M x 1 vector f of function evaluations
    d = dividedDiffNewton(x,y)
    N = len(x)
    M = len(z)
    f = zeros(M,1)
    for k in range(0,N):
        omega = ones(M,1)
        for j in range(0,k):
            for i in range(0,M):
                omega[i] = omega[i]*(z[i]-x[j])
        for i in range(0,M):
            f[i] = f[i] + omega[i]*d[k][k]
    return f

def interpLagrange2D(x,y,z,xval,yval):
    #  Input: N x 1 vector x of interpolation nodes
    #         M x 1 vector y of interpolation nodes
    #         N x M vector z of values at interpolation nodes
    #         P x 1 vector xval of nodes for function evaluation
    #         R x 1 vector yval of nodes for function evaluation
    #  Output: P x R vector f of function evaluations
    N = len(x)
    M = len(y)
    P = len(xval)
    R = len(yval)
    f = zeros(P,R)
    for r in range(0,R):
        for p in range(0,P):
            fpr = 0.0
            for i in range(0,N):
                for j in range(0,M):
                    fpr = fpr + z[i][j]*interpolantLagrange(x,xval[p],i)[0]*interpolantLagrange(y,yval[r],j)[0]
            f[p][r] = fpr
    return f

def interpLagrange3D(x,y,z,w,xval,yval,zval):
    #  Input: N x 1 vector x of interpolation nodes
    #         M x 1 vector y of interpolation nodes
    #         L x 1 vector z of interpolation nodes
    #         L x M x N vector w of values at interpolation nodes
    #         P x 1 vector xval of nodes for function evaluation
    #         R x 1 vector yval of nodes for function evaluation
    #         Q x 1 vector zval of nodes for function evaluation
    #  Output: Q x R x P vector f of function evaluations
    N = len(x)
    M = len(y)
    L = len(z)
    P = len(xval)
    R = len(yval)
    Q = len(zval)

    f = zeroTensor([Q,R,P])

    for q in range(0,Q):
        for r in range(0,R):
            for p in range(0,P):
                fprq = 0.0
                for k in range(0,L):
                    for j in range(0,M):
                        for i in range(0,N):
                            fprq = fprq + w[k][j][i]*interpolantLagrange(x,xval[p],i)[0]*interpolantLagrange(y,yval[r],j)[0]*interpolantLagrange(z,zval[q],k)[0]
                f[q][r][p] = fprq
    return f

def main(argv):



if __name__ == "__main__":
    main(sys.argv[1:])

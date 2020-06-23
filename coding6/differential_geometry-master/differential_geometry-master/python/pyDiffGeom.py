# -*- coding: utf-8 -*-

'''
=====================================================================================

Copyright (c) 2018 Université de Lorraine & Luleå tekniska universitet
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
from ..finite_differences/python import computeWeights, computeDerivativeAtPoint
from ..interpolation/python import zeros, zeroTensor

def det(A): # determinant of NxN matrix
    N = len(A)
    for row in A:
        if N!=len(row):
            raise TypeError('The matrix is not square!')
    res = 0.0
    if N==2:
        res = A[0][0]*A[1][1]-A[1][0]*A[0][1]
    else:
        for j in range(0,N):
            M = []
            for row in A(1:):
                row = []
                for v,value in enumerate(row):
                    if v!=j:
                        row.append(value)
                M.append(row)
            res += ((-1)**(2+j))*det(M)
    return res

def convertIndecesTensorToHelical(indeces,lengths):
    # indeces = [innermost ... outermost]
    # lengths = [innermost ... outermost]
    hindex = 0
    for i in range(len(indeces)-1,-1,-1):
        coeff = 1
        for j in range(i-1,-1,-1):
            coeff *= lengths[j]
        index += indeces[i]*coeff
    return hindex

def buildStencilsIndeces(indeces,lengths,stencilSize):
    # finite differences are assumed to be centered, thus stencilSize must an odd number
    # if the central finite difference cannot be performed, an alternative is computed
    stencilIndeces = []
    deltaIndex = np.floor(0.5*stencilSize)
    for i,cartesianIndex in enumerate(indeces):
        if (cartesianIndex-deltaIndex)>-1 && (cartesianIndex+deltaIndex)<lengths[i]: #centered finite difference can be performed
            startIndex = cartesianIndex-deltaIndex
            endIndex = cartesianIndex+deltaIndex
        elif (cartesianIndex-deltaIndex)>-1 # got to be ajusted on the right side
            startIndex = (cartesianIndex-deltaIndex) - ((cartesianIndex+deltaIndex)-(lengths[i]-1))
            endIndex = lengths[i]-1
            if startIndex<0:
                raise ValueError('The specified stencil size is not compatible with the mesh along dimension ' + str(i+1) + '. Check your input.')
        elif (cartesianIndex+deltaIndex)<lengths[i] # got to be ajusted on the left side
            startIndex = 0
            endIndex = (cartesianIndex+deltaIndex) - (cartesianIndex-deltaIndex) #(cartesianIndex-deltaIndex)<0 in this case
            if endIndex>=lengths[i]:
                raise ValueError('The specified stencil size is not compatible with the mesh along dimension ' + str(i+1) + '. Check your input.')
        else:
            raise ValueError('The specified stencil size is not compatible with the mesh along dimension ' + str(i+1) + '. Check your input.')
        currentStencil = []
        tempCartesianIndeces = indeces
        for j in range(startIndex,endIndex+1):
            tempCartesianIndeces[i] = j
            currentStencil.append(convertIndecesTensorToHelical(tempCartesianIndeces,lengths))
        stencilIndeces.append(currentStencil)
    return stencilIndeces

def covariantBaseAtPoint(indeces,lengths,rmap,qmap,stencilSize):
    index = convertIndecesTensorToHelical(indeces,lengths)
    rs = rmap[index]
    qs = qmap[index]
    hStencils = buildStencilsIndeces(indeces,lengths,stencilSize)
    gs = []
    for q,q0 in enumerate(qs):
        stencilIndexSet = hStencils[q]
        xs = []
        for sIndex in stencilIndexSet:
            xs.append(qmap[sIndex][q])
        g = []
        for r,rCoord in enumerate(rs):
            fs = []
            for sIndex in stencilIndexSet:
                fs.append(rmap[sIndex][r])
            g.append(computeDerivativeAtPoint(1,q0,xs,fs))
        gs.append(g)
    return gs

def contravariantBaseAtPoint(indeces,lengths,rmap,qmap,stencilSize):
    index = convertIndecesTensorToHelical(indeces,lengths)
    rs = rmap[index]
    qs = qmap[index]
    hStencils = buildStencilsIndeces(indeces,lengths,stencilSize)
    gs = []
    for r,r0 in enumerate(rs):
        stencilIndexSet = hStencils[r]
        xs = []
        for sIndex in stencilIndexSet:
            xs.append(rmap[sIndex][r])
        g = []
        for q,qCoord in enumerate(qs):
            fs = []
            for sIndex in stencilIndexSet:
                fs.append(qmap[sIndex][q])
            g.append(computeDerivativeAtPoint(1,r0,xs,fs))
        gs.append(g)
    return gs

def metricTensorAtPoint(gs):
    N = len(gs)
    M = len(gs[0])
    gij = zeroTensor([N,N])
    for i in range(0,N):
        for j in range(0,N):
            gValue = 0.0
            for k in range(0,M):
                gValue += gs[i][k]*gs[j][k]
            g[i][i] = gValue
    return gij

def computeVectorComponentsAtPoint(v,gs):
    N = len(gs)
    M = len(v)
    w = zeros(N,1)
    for i in range(0,N):
        value = 0.0
        for k in range(0,M):
            value += v[k]*gs[i][k]
        w[i] = value
    return w

def computeChristoffelSymbolAtPoint(indeces,lengths,qmap,covBase,gs,stencilSize):
    index = convertIndecesTensorToHelical(indeces,lengths)
    qs = qmap[index]
    hStencils = buildStencilsIndeces(indeces,lengths,stencilSize)
    K = len(gs)
    I = len(covBase[index])
    gammakij = zeroTensor([K,I,I])
    for k in range(0,K):
        for i in range(0,I):
            for j in range(0,I):
                q0 = qs[j]
                stencilIndexSet = hStencils[j]
                xs = []
                for sIndex in stencilIndexSet:
                    xs.append(qmap[sIndex][j])
                dg = []
                for r,rCoord in enumerate(covBase[index][i]):
                    fs = []
                    for sIndex in stencilIndexSet:
                        fs.append(covBase[sIndex][i][r])
                    dg.append(computeDerivativeAtPoint(1,q0,xs,fs))
                value = 0.0
                for l in range(0,I):
                    value += dg[l]*gs[k][l]
                gammakij[k][i][j] = value
    return gammakij

def computeContravariantRiemannAtPoint(indeces,lengths,qmap,connectionCoeffs,stencilSize):
    index = convertIndecesTensorToHelical(indeces,lengths)
    qs = qmap[index]
    hStencils = buildStencilsIndeces(indeces,lengths,stencilSize)
    L = len(connectionCoeffs)
    I = len(connectionCoeffs[0])
    J = I
    K = K
    S = len(connectionCoeffs[0][0])
    Rlijk = zeroTensor([L,I,J,K])
    for l in range(0,L):
        for i in range(0,I):
            for j in range(0,J):
                for k in range(0,K):
                    element = 0.0
                    # compute gamma^l_ik,j
                    dgamma = 0.0
                    q0 = qs[j]
                    stencilIndexSet = hStencils[j]
                    xs = []
                    for sIndex in stencilIndexSet:
                        xs.append(qmap[sIndex][j])
                    fs = []
                    for sIndex in stencilIndexSet:
                        fs.append(connectionCoeffs[sIndex][l][i][k])
                    dgamma = computeDerivativeAtPoint(1,q0,xs,fs)
                    element += dgamma
                    # compute gamma^l_ij,k
                    dgamma = 0.0
                    q0 = qs[k]
                    stencilIndexSet = hStencils[k]
                    xs = []
                    for sIndex in stencilIndexSet:
                        xs.append(qmap[sIndex][k])
                    fs = []
                    for sIndex in stencilIndexSet:
                        fs.append(connectionCoeffs[sIndex][l][i][j])
                    dgamma = computeDerivativeAtPoint(1,q0,xs,fs)
                    element -= dgamma
                    # compute gamma^l_js * gamma^s_ik
                    for s in range(0,S):
                        element += connectionCoeffs[index][l][j][s]*connectionCoeffs[index][s][i][k]
                    # compute gamma^l_ks * gamma^s_ij
                    for s in range(0,S):
                        element -= connectionCoeffs[index][l][k][s]*connectionCoeffs[index][s][i][j]
                    Rlijk[l][i][j][k] = element
    return Rlijk

def computeCovarRiemannAtPoint(gis,Rsklm):
    S = len(Rsklm)
    I = len(gis)
    K = len(Rsklm[0])
    L = len(Rsklm[0][0])
    M = len(Rsklm[0][0][0])
    Riklm = zeroTensor([I,K,L,M])
    for i in range(0,I):
        for k in range(0,K):
            for l in range(0,L):
                for m in range(0,M):
                    element = 0.0
                        for s in range(0,S):
                            element += gis[i][s]*Rsklm[s][k][l][m]
                    Riklm[i][k][l][m] = element
    return Riklm

def computeRicciTensorAtPoint(glm,Riljm):
    I = len(Riljm)
    L = len(Riljm[0])
    J = len(Riljm[0][0])
    M = len(Riljm[0][0][0])
    Rij = zeroTensor([I,J])
    for i in range(0,I):
        for j in range(0,J):
            element = 0.0
            for l in range(0,L):
                for m in range(0,M):
                    element += glm[l][m]*Riljm[i][l][j][m]
            Rij[i][j] = element
    return Rij

def computeRicciScalarAtPoint(gij,Rij):
    I = len(gij)
    J = len(gij[0])
    R = 0.0
    for i in range(0,I):
        for j in range(0,J):
            R += gij[i][j]*Rij[i][j]
    return R

def main(argv):



if __name__ == "__main__":
    main(sys.argv[1:])

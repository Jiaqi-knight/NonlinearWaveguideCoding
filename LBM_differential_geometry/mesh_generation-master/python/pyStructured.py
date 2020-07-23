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

def transfiniteLagrangeInterp2D(vertices,edges):
    # vertices = [c1,c2,c3,C4]
    # edges = [e1,e2,e3,e4]
    # e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
    # e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
    # e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
    # e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
    # c1 = c1(xi_min,eta_min) = (x1,y1)
    # c2 = c2(xi_max,eta_min) = (x2,y2)
    # c3 = c3(xi_max,eta_max) = (x3,y3)
    # c4 = c4(xi_min,eta_max) = (x4,y4)
    # e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    # e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    # e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    # e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    # c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    # c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    # c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    # c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    R = len(edges['e1'])
    P = len(edges['e2'])
    for r in range(0,R):
        for p in range(0,P):

def transfiniteLagrangeInterp3D():
    #                 f1 = e1(xi,eta_min,zeta) = f1(xi) = (x(xi),y(xi))
    #                 f2 = e2(xi_min,eta,zeta) = f2(eta) = (x(eta),y(eta))
    #                 f3 = e3(xi,eta_max,zeta) = f3(xi) = (x(xi),y(xi))
    #                 f4 = e4(xi_max,eta,zeta) = f4(eta) = (x(eta),y(eta))
    #                 f5 = e1(xi,eta_min,zeta) = f5(xi) = (x(xi),y(xi))
    #                 f6 = e2(xi_min,eta,zeta) = f6(eta) = (x(eta),y(eta))
    #                 e1 = e1(xi,eta_min,zeta) = e1(xi) = (x(xi),y(xi))
    #                 e2 = e2(xi_min,eta,zeta) = e2(eta) = (x(eta),y(eta))
    #                 e3 = e3(xi,eta_max,zeta) = e3(xi) = (x(xi),y(xi))
    #                 e4 = e4(xi_max,eta,zeta) = e4(eta) = (x(eta),y(eta))
    #                 e5 = e1(xi,eta_min,zeta) = e5(xi) = (x(xi),y(xi))
    #                 e6 = e2(xi_min,eta,zeta) = e6(eta) = (x(eta),y(eta))
    #                 e7 = e3(xi,eta_max,zeta) = e7(xi) = (x(xi),y(xi))
    #                 e8 = e4(xi_max,eta,zeta) = e8(eta) = (x(eta),y(eta))
    #                 e9 = e1(xi,eta_min,zeta) = e9(xi) = (x(xi),y(xi))
    #                 e10 = e2(xi_min,eta,zeta) = e10(eta) = (x(eta),y(eta))
    #                 e11 = e3(xi,eta_max,zeta) = e11(xi) = (x(xi),y(xi))
    #                 e12 = e4(xi_max,eta,zeta) = e12(eta) = (x(eta),y(eta))
    #                 c1 = c1(xi_min,eta_min) = (x1,y1)
    #                 c2 = c2(xi_max,eta_min) = (x2,y2)
    #                 c3 = c3(xi_max,eta_max) = (x3,y3)
    #                 c4 = c4(xi_min,eta_max) = (x4,y4)
    #                 c5 = c1(xi_min,eta_min) = (x1,y1)
    #                 c6 = c2(xi_max,eta_min) = (x2,y2)
    #                 c7 = c3(xi_max,eta_max) = (x3,y3)
    #                 c8 = c4(xi_min,eta_max) = (x4,y4)
    #                 f1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 f2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 f3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 f4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 f5 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 f6 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e5 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e6 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e7 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e8 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e9 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e10 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e11 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 e12 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
    #                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c5 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c6 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c7 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
    #                 c8 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)

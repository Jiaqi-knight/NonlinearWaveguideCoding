function [psi, gpsi, hs] = elm6_interp ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
   ,al,be,ga, xi,eta)
                                        
%==================================================
% Evaluation the surface metric coefficient, $h_s$,
% and computation of the basis functions and
% their gradients, over a 6-node triangle.
%==================================================

%--------
% prepare
%--------

alc = 1.0-al;
bec = 1.0-be;
gac = 1.0-ga;

alalc = al*alc;
bebec = be*bec;
gagac = ga*gac;

%----------------------------
% compute the basis functions
%----------------------------

psi(2) = xi*(xi-al+eta*(al-ga)/gac)/alc;
psi(3) = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec;
psi(4) = xi*(1.0-xi-eta)/alalc;
psi(5) = xi*eta/gagac;
psi(6) = eta*(1.0-xi-eta)/bebec;
psi(1) = 1.0-psi(2)-psi(3)-psi(4)-psi(5)-psi(6);
      
%--------------------------------------------------
% compute the xi derivatives of the basis functions
%--------------------------------------------------

dps2 =  (2.0*xi-al+eta*(al-ga)/gac)/alc;
dps3 =  eta*(be+ga-1.0)/(ga*bec);
dps4 =  (1.0-2.0*xi-eta)/alalc;
dps5 =  eta/gagac;
dps6 = -eta/bebec;
dps1 = -dps2-dps3-dps4-dps5-dps6;

%---------------------------------------------------
% compute the eta derivatives of the basis functions
%---------------------------------------------------

pps2 =  xi*(al-ga)/(alc*gac);
pps3 =  (2.0*eta-be+xi*(be+ga-1.0)/ga)/bec;
pps4 = -xi/alalc;
pps5 =  xi/gagac;
pps6 =  (1.0-xi-2.0*eta)/bebec;
pps1 = -pps2-pps3-pps4-pps5-pps6;

%----------------------------------------
% compute the xi and eta derivatives of x
%----------------------------------------

DxDxi = x1*dps1 + x2*dps2 + x3*dps3 + x4*dps4 + x5*dps5 + x6*dps6;
DyDxi = y1*dps1 + y2*dps2 + y3*dps3 + y4*dps4 + y5*dps5 + y6*dps6;

DxDet = x1*pps1 + x2*pps2 + x3*pps3 + x4*pps4 + x5*pps5 + x6*pps6;
DyDet = y1*pps1 + y2*pps2 + y3*pps3 + y4*pps4 + y5*pps5 + y6*pps6;

%------------------------------
% compute the surface metric hs
%------------------------------

vnz = DxDxi * DyDet - DxDet * DyDxi;
hs = sqrt(vnz^2);

%-----------------------------------
% compute the gradient of the basis functions
% by solving two linear equations:
%
% dx/dxi . grad = d psi/dxi
% dx/det . grad = d psi/det
%
% The system is solved by Cramer's rule
%-----------------------------------

A11 = DxDxi; A12 = DyDxi;
A21 = DxDet; A22 = DyDet;

Det = A11*A22-A21*A12;

B1 = dps1; B2 = pps1;
Det1 =   B1*A22 - B2*A12;
Det2 = - B1*A21 + B2*A11;
gpsi(1,1) = Det1/Det;
gpsi(1,2) = Det2/Det;

B1 = dps2; B2 = pps2;

Det1 =   B1*A22 - B2*A12;
Det2 = - B1*A21 + B2*A11;
gpsi(2,1) = Det1/Det;
gpsi(2,2) = Det2/Det;

B1 = dps3; B2 = pps3;

Det1 = B1*A22 - B2*A12;
Det2 = - B1*A21 + B2*A11;
gpsi(3,1) = Det1/Det;
gpsi(3,2) = Det2/Det;

B1 = dps4; B2 = pps4;

Det1 = B1*A22 - B2*A12; Det2 = - B1*A21 + B2*A11;
gpsi(4,1) = Det1/Det; gpsi(4,2) = Det2/Det;

B1 = dps5; B2 = pps5;

Det1 = B1*A22 - B2*A12; Det2 = - B1*A21 + B2*A11;
gpsi(5,1) = Det1/Det; gpsi(5,2) = Det2/Det;

B1 = dps6; B2 = pps6;

Det1 = B1*A22 - B2*A12; Det2 = - B1*A21 + B2*A11;
gpsi(6,1) = Det1/Det; gpsi(6,2) = Det2/Det;

%-----
% done
%-----

return;

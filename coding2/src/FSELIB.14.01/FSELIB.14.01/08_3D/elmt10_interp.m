function [psi, gpsi, hv] = elmt10_interp ...
...
   (x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4 ...
   ,x5,y5,z5, x6,y6,z6 ...
   ,x7,y7,z7, x8,y8,z8, x9,y9,z9, x10,y10,z10 ...
   ,xi,et,zt)

%==================================================
% Evaluation the surface metric coefficient, $h_s$,
% and computation of the basis functions and
% their gradients, over a 10-node tetrahedron
%==================================================

%--------
% prepare
%--------

om = 1.0-xi-et-zt;

%----------------------------
% compute the basis functions
%----------------------------

psi(1) = om*(2*om-1.0);
psi(2) = xi*(2*xi-1.0);
psi(3) = et*(2*et-1.0);
psi(4) = zt*(2*zt-1.0);
psi(5) = 4.0*xi*om;
psi(6) = 4.0*xi*et;
psi(7) = 4.0*et*om;
psi(8) = 4.0*zt*om;
psi(9) = 4.0*xi*zt;
psi(10)= 4.0*et*zt;

%--------------------------------------------------
% compute the xi derivatives of the basis functions
%--------------------------------------------------

dps(1) = -4*om+1.0;
dps(2) =  4*xi-1.0;
dps(3) =  0.0;
dps(4) =  0.0;
dps(5) =  4*om-4*xi;
dps(6) =  4*et;
dps(7) = -4*et;
dps(8) = -4*zt;
dps(9) =  4*zt;
dps(10)=  0.0;

%---------------------------------------------------
% compute the eta derivatives of the basis functions
%---------------------------------------------------

pps(1) = -4*om+1.0;
pps(2) = 0.0;
pps(3) = 4*et-1.0;
pps(4) = 0.0;
pps(5) = -4*xi;
pps(6) =  4*xi;
pps(7) = 4*om-4*et;
pps(8) = -4*zt;
pps(9) = 0.0;
pps(10)= 4*zt;

%----------------------------------------------------
% compute the zeta derivatives of the basis functions
%----------------------------------------------------

qps(1) = -4*om+1.0;
qps(2) = 0.0;
qps(3) = 0.0;
qps(4) = 4*zt-1.0;
qps(5) = -4*xi;
qps(6) = 0.0;
qps(7) = -4*et;
qps(8) = 4*om-4*zt;
qps(9) = 4*xi;
qps(10)= 4*et;

%sum1 = 0;
%sum2 = 0;
%sum3 = 0;
%sum4 = 0;
%for i=1:10
% sum1 = sum1+psi(i);
% sum2 = sum2+dps(i);
% sum3 = sum3+pps(i);
% sum4 = sum4+qps(i);
%end
%[sum1 sum2 sum3 sum4]
%pause

%--------------------------------
% compute the xi derivatives of x
%--------------------------------

DxDxi = x1*dps(1) + x2*dps(2) + x3*dps(3) ...
      + x4*dps(4) + x5*dps(5) + x6*dps(6) ...
      + x7*dps(7) + x8*dps(8) + x9*dps(9) + x10*dps(10);

DyDxi = y1*dps(1) + y2*dps(2) + y3*dps(3) ...
      + y4*dps(4) + y5*dps(5) + y6*dps(6) ...
      + y7*dps(7) + y8*dps(8) + y9*dps(9) + y10*dps(10);

DzDxi = z1*dps(1) + z2*dps(2) + z3*dps(3) ...
      + z4*dps(4) + z5*dps(5) + z6*dps(6) ...
      + z7*dps(7) + z8*dps(8) + z9*dps(9) + z10*dps(10);

%--------------------------------
% compute the eta derivatives of x
%--------------------------------

DxDet = x1*pps(1) + x2*pps(2) + x3*pps(3) ...
      + x4*pps(4) + x5*pps(5) + x6*pps(6) ...
      + x7*pps(7) + x8*pps(8) + x9*pps(9) + x10*pps(10);

DyDet = y1*pps(1) + y2*pps(2) + y3*pps(3) ...
      + y4*pps(4) + y5*pps(5) + y6*pps(6) ...
      + y7*pps(7) + y8*pps(8) + y9*pps(9) + y10*pps(10);

DzDet = z1*pps(1) + z2*pps(2) + z3*pps(3) ...
      + z4*pps(4) + z5*pps(5) + z6*pps(6) ...
      + z7*pps(7) + z8*pps(8) + z9*pps(9) + z10*pps(10);

%--------------------------------
% compute the zeta derivatives of x
%--------------------------------

DxDzt = x1*qps(1) + x2*qps(2) + x3*qps(3) ...
      + x4*qps(4) + x5*qps(5) + x6*qps(6) ...
      + x7*qps(7) + x8*qps(8) + x9*qps(9) + x10*qps(10);

DyDzt = y1*qps(1) + y2*qps(2) + y3*qps(3) ...
      + y4*qps(4) + y5*qps(5) + y6*qps(6) ...
      + y7*qps(7) + y8*qps(8) + y9*qps(9) + y10*qps(10);

DzDzt = z1*qps(1) + z2*qps(2) + z3*qps(3) ...
      + z4*qps(4) + z5*qps(5) + z6*qps(6) ...
      + z7*qps(7) + z8*qps(8) + z9*qps(9) + z10*qps(10);

%------------------------------
% compute the surface metric hs
%------------------------------

Jac = [DxDxi DxDet DxDzt;
       DyDxi DyDet DyDzt;
       DzDxi DzDet DzDzt];

hv = det(Jac);

%-----------------------------------
% compute the gradient of the basis functions
% by solving two linear equations:
%
% dx/dxi . grad = d psi/dxi
% dx/det . grad = d psi/deta
% dx/dzt . grad = d psi/dzeta
%-----------------------------------

 for i=1:10
   rhs = [dps(i), pps(i), qps(i)];
   sln = rhs/Jac;
   gpsi(i,1) = sln(1);
   gpsi(i,2) = sln(2);
   gpsi(i,3) = sln(3);
 end

%-----
% done
%-----

return;

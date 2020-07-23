function [ap,bp,cp,dp,ep,b] = sdq_sys (ne,xe,q0,fL,k,s)

%======================================================
% Compact assembly of the pentadiagonal linear system 
% for one-dimensional steady diffusion
% with quadratic elements (sdq)
%======================================================

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%------------------------------
% number of unique global nodes
%------------------------------

ng = 2*ne+1;

%------------------------------------
% initialize the pentadiagonal matrix
%------------------------------------

ap = zeros(ng,1); bp = zeros(ng,1); cp = zeros(ng,1);
dp = zeros(ng,1); ep = zeros(ng,1);

%-------------------------------
% initialize the right-hand side
%-------------------------------

b = zeros(ng,1); b(1) = q0/k;

%-----------------------
% loop over all elements
%-----------------------

for l=1:ne

  cf = 1.0/(3.0*h(l));
  A11 = 7.0*cf;  A12 = -8.0*cf; A13 = cf;
  A21 = A12;     A22 = 16.0*cf; A23 = A12;
  A31 = A13;     A32 = A23;     A33 = A11;

  cf = h(l)/30.0;
  B11 = 4.0*cf; B12 =  2.0*cf;  B13 = -cf;
  B21 = B12;    B22 = 16.0*cf;  B23 = B12;
  B31 = B13;    B32 = B23;      B33 = B11;

  cl1 = 2*l-1; cl2 = 2*l; cl3 = 2*l+1;

  ap(cl1) = ap(cl1) + A11;
  bp(cl1) = bp(cl1) + A12;  
  cp(cl1) = cp(cl1) + A13;  

  dp(cl2) = dp(cl2) + A21;
  ap(cl2) = ap(cl2) + A22;
  bp(cl2) = bp(cl2) + A23;

  ep(cl3) = ep(cl3) + A31;
  dp(cl3) = dp(cl3) + A32;
  ap(cl3) = ap(cl3) + A33;

  b(cl1) = b(cl1) + (B11*s(cl1) + B12*s(cl2) + B13*s(cl3))/k;
  b(cl2) = b(cl2) + (B21*s(cl1) + B22*s(cl2) + B23*s(cl3))/k;
  b(cl3) = b(cl3) + (B31*s(cl1) + B32*s(cl2) + B33*s(cl3))/k;

end

%----------------------------------
% implement the Dirichlet condition
%----------------------------------

b(ng-2) = b(ng-2) - cp(ng-2)*fL;
b(ng-1) = b(ng-1) - bp(ng-1)*fL;

%-----
% done
%-----

return;

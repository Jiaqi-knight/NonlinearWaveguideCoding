function [at,bt,ct,b] ...
...
   = cdl_sys (ne,xe,q0,ht,finf,fL,k,U,rho,cp,s,leftbc)

%--------------------------------------------------
% Compact assembly of the tridiagonal linear system
% for steady one-dimensional convection--diffusion
% with linear elements (scdl)
%--------------------------------------------------

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%-----------
% initialize
%-----------

at = zeros(ne,1);
bt = zeros(ne,1);
ct = zeros(ne,1);
b = zeros(ne,1);

if(leftbc==1)       % Neumann
 b(1) = q0/k;
elseif(leftbc==2)   % Robin
 b(1) = ht*finf/k;
 at(1) = ht/k;
end

cf = 0.5*U*rho*cp/k;

C11 = -cf; C12 = cf;    % advection matrix
C21 = -cf; C22 = cf;

%------------------------
% loop over ne-1 elements
%------------------------

for l=1:ne-1

  A11 = 1/h(l); A12=-A11; A22=A11;
  B11 = h(l)/3.0; B12=0.5*B11; B22=B11;

  at(l)   = at(l)   + A11 + C11;
  bt(l)   = bt(l)   + A12 + C12;  
  ct(l+1) = ct(l+1) + A12 + C21; 
  at(l+1) = at(l+1) + A22 + C22;
  b(l)    = b(l)   + (B11*s(l) + B12*s(l+1))/k;
  b(l+1)  = b(l+1) + (B12*s(l) + B22*s(l+1))/k;

end

%-------------
% last element
%-------------

A11 = 1/h(ne); A12=-A11;
B11 = h(ne)/3.0; B12=0.5*B11;

at(ne) = at(ne) + A11+C11;
b(ne) = b(ne) + (B11*s(ne) + B12*s(ne+1)) /k - (A12+C12)*fL;

%-----
% done
%-----

return;

function [at,bt,ct,b] = hlml_sys (ne,xe,q0,fL,k,alpha)

%==================================================
% FSELIB
%
% Assembly of the tridiagonal linear system 
% for the one-dimensional Helmholtz equation 
%==================================================

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%----------------------------------
% initialize the tridiagonal matrix
%----------------------------------

at = zeros(ne,1);
bt = zeros(ne,1);
ct = zeros(ne,1);

%-------------------------------
% initialize the right-hand side
%-------------------------------

b = zeros(ne,1);

b(1) = q0/k;

%----------------------------------
% loop over the first ne-1 elements
%----------------------------------

for l=1:ne-1

  A11 = 1/h(l); A12 =-A11;
  A21 = A12;    A22 = A11;

  B11 = h(l)/3.0; B12 = 0.5*B11;
  B21 = B12;      B22 = B11;

  at(l) = at(l) + A11 - alpha  * B11;
  bt(l) = bt(l) + A12 - alpha * B12;  

  ct(l+1) = ct(l+1) + A21 - alpha * B21; 
  at(l+1) = at(l+1) + A22 - alpha * B22;

end

%----------------------------
% the last element is special
%----------------------------

A11 = 1.0/h(ne); A12 =-A11;

B11 = h(ne)/3.0; B12 = 0.5*B11;

at(ne) = at(ne) + A11 - alpha * B11;
b(ne) = b(ne) - (A12- alpha * B12)*fL;

%-----
% done
%-----

return;

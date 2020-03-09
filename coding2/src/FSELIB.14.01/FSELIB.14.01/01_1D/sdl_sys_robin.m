function [at,bt,ct,b] = sdl_sys_robin (ne,xe,q0,fL,k,s,finfty,hT)

%==================================================
% FSELIB
%
% compact assembly of a tridiagonal linear system 
% for one-dimensional steady diffusion
% with linear elements (sdl)
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

b(1) = hT*finfty/k
at(1) = hT/k

%----------------------------------
% loop over the first ne-1 elements
%----------------------------------

for l=1:ne-1

  A11 = 1/h(l); A12 =-A11;
  A21 = A12;    A22 = A11;

  B11 = h(l)/3.0; B12 = 0.5*B11;
  B21 = B12;      B22 = B11;

  at(l) = at(l) + A11;
  bt(l) = bt(l) + A12;  

  ct(l+1) = ct(l+1) + A21; 
  at(l+1) = at(l+1) + A22;

  b(l)   = b(l)   + (B11*s(l) + B12*s(l+1))/k;
  b(l+1) = b(l+1) + (B21*s(l) + B22*s(l+1))/k;

end

%----------------------------
% the last element is special
%----------------------------

A11 = 1.0/h(ne); A12 =-A11;

B11 = h(ne)/3.0; B12 = 0.5*B11;

at(ne) = at(ne) + A11;
b(ne) = b(ne) + (B11*s(ne) + B12*s(ne+1)) /k - A12*fL;

%-----
% done
%-----

return;

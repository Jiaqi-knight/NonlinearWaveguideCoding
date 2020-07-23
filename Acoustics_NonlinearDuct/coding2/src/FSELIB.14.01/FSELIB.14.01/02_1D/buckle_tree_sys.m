function [at,bt,ct,ar,br,cr] = buckle_tree_sys (ne,xe,L)

%==================================================
% FSELIB
%
% compact assembly of a tridiagonal linear system 
% for the tree buckling equation
%
% at,bt,ct correspond to the diffusion matix
% ar,br,cr correspond to the generalized mass matix
%==================================================

%---------------------------
% element size and midpoints
%---------------------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
  xm(l) = 0.5*(xe(l+1)+xe(l));
end

%------------------------------------
% initialize the tridiagonal matrices
%------------------------------------

at = zeros(ne,1);
bt = zeros(ne,1);
ct = zeros(ne,1);

ar = zeros(ne,1);
br = zeros(ne,1);
cr = zeros(ne,1);

%-------------------------------
% the first element is special
% due to the Dirichlet condition
%-------------------------------

A11 = 1.0/h(1); A12 =-A11;

B11 = h(1)/3.0; B12 = 0.5*B11;

at(1) = at(1) + A11;
xx = L-xm(1);
ar(1) = ar(1) + xx * B11;

%---------------------------------------
% loop over the subsequent ne-1 elements
%---------------------------------------

for l=2:ne

  A11 = 1/h(l); A12 =-A11;
  A21 = A12;    A22 = A11;

  at(l-1) = at(l-1) + A11;
  bt(l-1) = bt(l-1) + A12;
  ct(l) = ct(l) + A21;
  at(l) = at(l) + A22;

  B11 = h(l)/3.0; B12 = 0.5*B11;
  B21 = B12;      B22 = B11;

  xx = L-xm(l);
  ar(l-1) = ar(l-1) + xx * B11;
  br(l-1) = br(l-1) + xx * B12;  
  cr(l) = cr(l) + xx * B21; 
  ar(l) = ar(l) + xx * B22;

end

%-----
% done
%-----

return;

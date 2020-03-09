function [at,bt,ct,b] = sdq_cnd_sys (ne,xe,q0,fL,k,s)

%============================================================
% Compact assembly of the condensed tridiagonal linear system 
% for one-dimensional steady diffusion
% with quadratic elements (sdqc)
%============================================================

%-------------
% element size
%-------------

for l=1:ne
  h(l)=xe(l+1)-xe(l);
end

%------------------------------------
% initialize the tridiagonal matrix
%------------------------------------

at = zeros(ne+1,1);
bt = zeros(ne+1,1);
ct = zeros(ne+1,1);

%-------------------------------
% initialize the right-hand side
%-------------------------------

b = zeros(ne+1,1);
b(1) = q0/k;

%-----------------------
% loop over all elements
%-----------------------

for l=1:ne

  cf = 1.0/(6.0*h(l));
  A11 = 14.0*cf; A12 =-16.0*cf; A13 = 2.0*cf;
  A21 = A12;     A22 = 32.0*cf; A23 = A12;
  A31 = A13;     A32 = A23;     A33 = A11;

  cf = h(l)/30.0;
  B11 = 4.0*cf; B12 =  2.0*cf;  B13 = -1.0*cf;
  B21 = B12;    B22 = 16.0*cf;  B23 = B12;
  B31 = B13;    B32 = B23;      B33 = B11;

  % condense the diffusion matrix:

  r12 = A12/A22;
  A11 = A11-r12*A21; A13= A13-r12*A23;

  r32 = A32/A22;
  A31 = A31-r32*A21; A33= A33-r32*A23;

  at(l)   = at(l)   + A11;
  bt(l)   = bt(l)   + A13;  
  ct(l+1) = ct(l+1) + A31;
  at(l+1) = at(l+1) + A33;

  % condense the right-hand side

  cl1=2*l-1; cl2=2*l; cl3=2*l+1;

  b(l)   = b(l)   +     (B11*s(cl1) + B12*s(cl2) + B13*s(cl3))/k;
  b(l)   = b(l)   - r12*(B21*s(cl1) + B22*s(cl2) + B23*s(cl3))/k;
  b(l+1) = b(l+1) +     (B31*s(cl1) + B32*s(cl2) + B33*s(cl3))/k;
  b(l+1) = b(l+1) - r32*(B21*s(cl1) + B22*s(cl2) + B23*s(cl3))/k;

end

%----------------------------------
% implement the Dirichlet condition
%----------------------------------

b(ne) = b(ne)-A13*fL;

%-----
% done
%-----

return;

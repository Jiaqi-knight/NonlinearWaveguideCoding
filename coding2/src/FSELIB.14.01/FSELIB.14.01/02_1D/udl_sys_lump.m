function [at,bt,ct,rhs] = udl_sys_lump (ne,xe,q0,f,fL,k,kappa,s,Dt)

%==================================================
% FSELIB
%
% Compact assembly of the tridiagonal linear system 
% for one-dimensional unsteady diffusion
% with linear elements (udl) and mass lumping
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

rhs = zeros(ne,1);
rhs(1) = Dt*q0/k;

%----------------------------------
% loop over the first ne-1 elements
%----------------------------------

for l=1:ne-1

  A11 = 1/h(l); A12 =-A11;      % diffusion matrix
  A21 = A12;    A22 = A11;

  B11 = h(l)/3.0; B12 = 0.5*B11;  % mass matrix
  B21 = B12;      B22 = B11;

  at(l)   = at(l) + B11;        % tridiagonal matrix components
  bt(l)   = bt(l) + B12;  
  ct(l+1) = ct(l+1) + B21; 
  at(l+1) = at(l+1) + B22;

%  rhs(l)   = rhs(l)   + B11*f(l) + B12*f(l+1);
%  rhs(l+1) = rhs(l+1) + B21*f(l) + B22*f(l+1);

  rhs(l)   = rhs(l)   + B11*f(l) + B12*f(l);
  rhs(l+1) = rhs(l+1) + B21*f(l+1) + B22*f(l+1);

  rhs(l)   = rhs(l)   - Dt*kappa*(A11*f(l) + A12*f(l+1));
  rhs(l+1) = rhs(l+1) - Dt*kappa*(A21*f(l) + A22*f(l+1));

  rhs(l)   = rhs(l)   + Dt*kappa*(B11*s(l) + B12*s(l+1))/k;
  rhs(l+1) = rhs(l+1) + Dt*kappa*(B12*s(l) + B22*s(l+1))/k;

end

%------------------------
% the last element is special
%------------------------

A11 = 1.0/h(ne); A12 =-A11;
B11 = h(ne)/3.0; B12 = 0.5*B11;

at(ne) = at(ne) + B11;

rhs(ne) = rhs(ne) + B11*f(ne);

rhs(ne) = rhs(ne) - Dt*kappa*(A11*f(ne) + A12*fL);
rhs(ne) = rhs(ne) + Dt*kappa*(B11*s(ne) + B12*s(ne+1))/k;

%-------------
% mass lumping
%-------------

at(1) = at(1)+bt(1);
bt(1) = 0;

for i=2:ne-1
 at(i) = at(i)+bt(i)+ct(i);
 bt(i) = 0;
 ct(i) = 0;
end

at(ne) = at(ne)+ct(ne);
ct(ne) = 0.0;

%-----
% done
%-----

return;

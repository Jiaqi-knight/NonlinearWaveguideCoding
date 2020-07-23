function [gdm,b] = sds_lob_sys (ne,xe,np,ng,c,q0,fL,k,s)

%=====================================
% Assembly of a linear system 
% for one-dimensional Steady Diffusion
% with Spectral elements (sds)
%
% gdm: global diffusion matrix
% b:   right-hand side
% c:   connectivity matrix
%=====================================

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%-----------
% initialize
%-----------

gdm = zeros(ng,ng);
b = zeros(1,ng);

b(1) = q0/k;    % Neumann condition on the left

%-----------------------
% loop over the elements
%-----------------------

for l=1:ne

   m = np(l);

   elm_mm = 0.5*h(l)*emm_lob_tbl(m);         % element mass matrix
%  elm_mm = 0.5*h(l)*emm_lob_lump_tbl(m);    % element mass matrix (lumped)
   elm_dm = 2.0*edm_lob_tbl(m)/h(l);         % element diffusion matrix

   for i=1:m+1
     i1 = c(l,i);
     for j=1:m+1
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + elm_dm(i,j);
       b(i1)      =   b(i1)     + elm_mm(i,j)*s(j1)/k;
     end
   end

end

%----------
% last node
%----------

m = np(ne);

for i=1:m
  i1 = c(ne,i);
  b(i1) = b(i1) - elm_dm(i,m+1)*fL;
end

%-----
% done
%-----

return;

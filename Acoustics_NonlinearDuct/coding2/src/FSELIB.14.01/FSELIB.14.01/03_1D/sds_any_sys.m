function [gdm,b] = sds_any_sys (ne,xe,np,ng,c,xien,q0,fL,k,s)

%==============================================
% Assemble a linear system
% for one-dimensional steady diffusion
% for any nodal distribution
%
% converted into the evenly-spaced nodal base
%
% gdm: global diffusion matrix
% b:   right-hand side
% c:   connectivity matrix
%==============================================

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
b   = zeros(1,ng);    % RHS

b(1) = q0/k;    % Neumann condition on the left

%-----------------------
% loop over the elements
%-----------------------

for l=1:ne

   m = np(l);

   clear xi

   for j=1:m+1
     xi(j) = xien(l,j);
   end

  vdm = vdm_modal_lob(m,xi); 
  elm_dm_modal = edm_modal_lob(m,h(l));
  elm_mm_modal = emm_modal_lob(m,h(l));
  elm_dm = inv(vdm)*elm_dm_modal*inv(vdm');
  elm_mm = inv(vdm)*elm_mm_modal*inv(vdm');

% assemble:

   for i=1:m+1
      i1 = c(l,i);
      for j=1:m+1
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + elm_dm(i,j);
       b(i1)      =   b(i1)    + elm_mm(i,j)*s(j1)/k;
      end
   end

%--
end   % end of loop over elements
%--

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

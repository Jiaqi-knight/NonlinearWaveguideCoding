function [gdm,b] = sds_modal_sys (ne,xe,np,ng,c,q0,fL,k,dg)

%=======================================
% assembly of a linear system
% for one-dimensional steady diffusion
% using the modal expansion
%=======================================

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%-----------
% initialize
%-----------

gdm = zeros(ng,ng); % global diffusion matrix
b = zeros(1,ng);    % RHS

b(1) = q0/k;    % Neumann condition on the left

%-----------------------
% loop over the elements
%-----------------------

for l=1:ne

   m = np(l);

% element diffusion and mass matrices

  clear elm_dm elm_mm
  elm_dm =  edm_modal_lob(m,h(l));
  elm_mm =  emm_modal_lob(m,h(l));

   for i=1:m+1
      i1 = c(l,i);
      for j=1:m+1
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + elm_dm(i,j);
       b(i1)      = b(i1)      + elm_mm(i,j)*dg(j1)/k;
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

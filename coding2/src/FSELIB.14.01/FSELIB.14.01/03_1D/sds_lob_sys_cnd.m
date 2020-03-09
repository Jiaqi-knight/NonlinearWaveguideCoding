function [gdm,b] = sds_lob_sys_cnd (ne,xe,np,c,q0,fL,k,s)

%--------------------------------------------------------------
% Assembly of the condensed linear system for
% one-dimensional Steady Diffusion with Spectral elements (sdsc)
%
% gdm: global diffusion matrix
% b:   right-hand side
% c:   connectivity matrix
%--------------------------------------------------------------

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%-----------
% initialize
%-----------

gdm = zeros(ne+1,ne+1);
b = zeros(1,ne+1);

b(1) = q0/k;

%-----------------------
% loop over the elements
%-----------------------

for l=1:ne

   m = np(l);   % element polynomial expansion order

   elm_dm = 2.0*edm_lob_tbl(m)/h(l);         % element diffusion matrix
   elm_mm = 0.5*h(l)*emm_lob_tbl(m);         % element mass matrix
%  elm_mm = 0.5*h(l)*emm_lob_lump_tbl(m);    % element mass matrix, lumped

   % Compute the element equations right-hand side

   for i=1:m+1
     belm(i) = 0.0;
     for j=1:m+1
      j1 = c(l,j);
      belm(i) = belm(i) + elm_mm(i,j)*s(j1);
     end
     belm(i) = belm(i)/k;
   end

   % Condense the element diffusion matrix
   % and right-hand side

%----------
   if(m > 1)
%----------

   for i=2:m      % loop over interior nodes
    for j=1:m+1   % loop over all element nodes

     if(i ~= j)   % skip the self node
       rji = elm_dm(j,i)/elm_dm(i,i);
       for p=1:m+1
        elm_dm(j,p) = elm_dm(j,p) - rji*elm_dm(i,p);
       end
       belm(j) = belm(j) - rji*belm(i);
     end

    end
   end

%----------
   end
%----------

% elm_dm  % print

%---
%  assign element equations to global equations
%---

   gdm(l,  l)   = gdm(l,l)      + elm_dm(1,  1);
   gdm(l,  l+1) = gdm(l, l+1)   + elm_dm(1,  m+1);
   gdm(l+1,l)   = gdm(l+1,l)    + elm_dm(m+1,1);
   gdm(l+1,l+1) = gdm(l+1, l+1) + elm_dm(m+1,m+1);

   b(l)   = b(l)   + belm(1);
   b(l+1) = b(l+1) + belm(m+1);

end

%----------
% last node
%----------

b(ne) = b(ne) - elm_dm(1,m+1)*fL;

%-----
% done
%-----

return;

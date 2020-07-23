function [gdm,gmm] = gdmm_lob (ne,xe,np,ng,c)

%=============================================
% assembly of the global diffusion matrix (gdm)
% and gobal mass matrix (gmm)
%
% c:  connectivity matrix
%=============================================

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
gmm = zeros(ng,ng);

%-----------------------
% loop over the elements
%-----------------------

for l=1:ne

   m = np(l);

    elm_mm = 0.5*h(l)*emm_lob_tbl(m);         % element mass matrix
%   elm_mm = 0.5*h(l)*emm_lob_lump_tbl(m);    % element mass matrix
    elm_dm = 2.0*edm_lob_tbl(m)/h(l);         % element diffusion matrix

   for ip=1:m+1
      i1 = c(l,ip);
      for jp=1:m+1
       i2 = c(l,jp);
       gdm(i1,i2) = gdm(i1,i2) + elm_dm(ip,jp);
       gmm(i1,i2) = gmm(i1,i2) + elm_mm(ip,jp);
      end
   end

end

%---
% done
%---

return;

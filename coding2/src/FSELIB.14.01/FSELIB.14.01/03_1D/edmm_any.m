clear all

%================================
% Code edmm_any
%
% element diffusion and mass
% matrices for arbitrary nodal
% discretization
%===============================

m = 4;
h = 1.0;

ichoose=1; % evenly spaced
ichoose=2; % lobatto

%---
% evenly spaced
%---

if(ichoose==1)
 for i=1:m+1
  xi(i) = -1+2*(i-1)/m;
 end
end

%---
% lobatto
%---

if(ichoose==2)
 xi(1) = -1.0;
 if(m>1)
   [zL, wL]= lobatto(m-1);
   for i=2:m
    xi(i) = zL(i-1);
   end
 end
 xi(m+1)=1.0;
end

%---
% compute the vandermond matrix
%---

vdm = vdm_modal_lob(m,xi);

%---
% compute the modal matrices
%---

elm_dm_modal = edm_modal_lob(m,h);
elm_mm_modal = emm_modal_lob(m,h);

%---
% compute the nodal matrices
%---

elm_dm = inv(vdm)*elm_dm_modal*inv(vdm');
elm_mm = inv(vdm)*elm_mm_modal*inv(vdm');

elm_dm
elm_mm

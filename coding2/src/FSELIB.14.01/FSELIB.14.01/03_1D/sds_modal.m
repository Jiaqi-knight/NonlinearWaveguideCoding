close all
clear all

%=======================================
% Code sds_modal
%
% modal expansion with spectral elements
%=======================================

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)
box on

%-----------
% input data
%-----------

L=1.0; k=1.0; fL=0; q0=-1.0;

ne=2; ratio=1.0;      % number of elements, stretch ratio
np(1) = 3;            % element order expansions
np(2) = 3;

%-----
% trap
%-----

for i=1:ne

   if(np(i)>7)
   disp(' -->');
   disp(' max polynomial order is 7');
   error(' sdms: Sorry this high order not yet implemented');
   end

end

%------------------------
% element node generation 
%------------------------

[xe,xen,xien,xg,c,ng] = discr_lob (0,L,ne,ratio,np);

%-------------------
% compute the source coefficients
%-------------------

for l=1:ne

  m = np(l);

  for j=1:m+1
   s(j) = 10.0*exp(-5.0*xen(l,j)^2/L^2);
  end

  vdm = vdm_modal_lob(m,xien(l,:));

  d = s/vdm;

  for j=1:m+1
    dg(c(l,j)) = d(j);
  end

  clear s vdm d

end

%-----------------
% element assembly
%-----------------

[gdm,b] = sds_modal_sys (ne,xe,np,ng,c,q0,fL,k,dg);

%--------------
% linear solver
%--------------

gdm(:,ng) = [];  % remove the last (ng) column
gdm(ng,:) = [];  % remove the last (ng) row
b(:,ng) = [];    % remove the last (ng) element

cg = b/gdm';       % solve the linear system
cg = [cg fL];      % add the value at the right end

%---
% extract the element modal coefficients
%---

shift = 0;
for l=1:ne
 for j=1:np(l)+1
  c(l,j) = cg(shift+j);
 end
 shift = shift+np(l);
end

%---
% reconstract the nodal values
%---

Ic=0;
for l=1:ne
  vdm = vdm_modal_lob(m,xien(l,:));
  felm = vdm'*c(l,:)';
  for j=1:np(l)
   Ic=Ic+1;
   f(Ic) = felm(j);
  end
end

f = [f fL];      % add the value at the right end

%-----
% plot
%-----

plot(xg, f,'k-+')

%-----
% done
%-----

clear all
close all

%========================================================
% Code sds_lob_cnd
%
% Steady one-dimensional diffusion with spectral elements
% using the condensed formulation where only the end-nodes
% appear in the final system of equations
%========================================================

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

L=1.0; k=1.0; q0=-1.0; fL=0.0;

ne=3; ratio=5.0;             % number of elements, stretch ratio
np(1)=5; np(2)=5; np(3)=3;   % element expansion order

%------------------------
% element node generation 
%------------------------

[xe,xen,xien,xg,c,ng] = discr_lob(0,L,ne,ratio,np);

%-------------------
% specify the source
%-------------------

for i=1:ng
 s(i) = 10.0*exp(-5.0*xg(i)^2/L^2);
end

%-----------------
% element assembly
%-----------------

[gdm,b] = sds_lob_sys_cnd (ne,xe,np,c,q0,fL,k,s);

%--------------
% linear solver
%--------------

gdm(:,ne+1) = [];  % remove the last (ne+1) column
gdm(ne+1,:) = [];  % remove the last (ne+1) row
b(:,ne+1) = [];    % remove the last (ne+1) node

f = b/gdm';       % solve the linear system

f = [f fL];      % add the value at the right end

%-----
% plot
%-----

plot(xe, f,'k-+')

ye = zeros(ne+1,1);
plot(xe,ye,'ko');
for i=1:ne
 plot([xe(i),xe(i)],[0,2.5],'k:');
end

%-----
% done
%-----

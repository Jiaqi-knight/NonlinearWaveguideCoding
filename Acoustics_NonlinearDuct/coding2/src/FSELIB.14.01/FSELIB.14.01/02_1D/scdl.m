clear all
%close all

%=====================================================
% FSELIB
%
% CODE scdl
%
% Code for steady one-dimensional convection-diffusion
% with linear elements (scdl)
%=====================================================

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)
box on

%-----------
% input data
%-----------

L    = 1.0;
k    = 1.0;
fL   = 0.0;
q0   =-1.0; 
ht   = 1.0;
finf = 0.0;

rho =1.0; cp=1.0;

ne = 16; ratio=1.0;

U =   0.0;
U =  20.0;
U =  60.0;
U =  80.0;
U = 100.0;

leftbc = 1; % newmann
leftbc = 2; % robin

%---------------
% grid generation
%---------------

xe=elm_line1 (0,L,ne,ratio);

%-------------------
% specify the source
%-------------------

for i=1:ne+1
 s(i) = 10.0*exp(-5.0*xe(i)^2/L^2);
end

%-----------------
% element assembly
%-----------------

[at,bt,ct,b] = scdl_sys (ne,xe,q0,ht,finf,fL,k,U,rho,cp,s,leftbc);

%--------------
% linear solver
%--------------

f = thomas(ne,at,bt,ct,b);

f(ne+1) = fL;

%-----
% plot
%-----

plot(xe, f,'k-o');

%-----
% done
%-----

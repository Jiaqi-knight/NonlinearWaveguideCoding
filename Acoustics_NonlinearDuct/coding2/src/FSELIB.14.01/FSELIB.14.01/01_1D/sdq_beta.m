clear all
close all

%=================================
% FSELIB
%
% Code sdq_beta
%
% steady one-dimensional diffusion
% with quadratic elements
% and arbitrary position of the
% interior node determined
% by the scalar parameter beta
%
% -1 < beta < 1
%
% beta=0 yields a midpoint node
%
% ne : number of elements
%=================================

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

L = 1.0; k = 1.0; q0 = -1.0; fL=0.0;
ne = 8; ratio = 5.0;
beta = 0.5;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%------------------------------
% number of unique global nodes
%------------------------------

ng = 2*ne+1;

%--------------------------
% generate the global nodes
%--------------------------

Ic = 0;  % counter

for i=1:ne
 Ic = Ic+1;
 xg(Ic) = xe(i);
 Ic = Ic+1;
 xg(Ic) = 0.5*(xe(i)+xe(i+1)) + 0.5*(xe(i+1)-xe(i))*beta;
end

xg(ng) = xe(ne+1);

%-------------------
% specify the source
%-------------------

for i=1:ng
 s(i) = 1.0;
 s(i) = 10.0*exp(-5.0*xg(i)^2/L^2);
end

%-------------------------
% compact element assembly
%-------------------------

[ap,bp,cp,dp,ep,b] = sdq_beta_sys (ne,xe,beta,q0,fL,k,s);

%--------------
% linear solver
%--------------

f = penta (ng-1,ap,bp,cp,dp,ep,b);

f(ng) = fL;

%-----
% plot
%-----

plot(xg, f,'-k+');

for i=1:ne+1
 plot([xe(i),xe(i)], [0, 2.5],'k:');
end

%-----
% done
%-----

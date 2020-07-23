clear all
close all

%============================================
% FSELIB
%
% Code sdq_cnd
%
% Steady one-dimensional diffusion
% with quadratic elements
% and a condensed formulation
%============================================

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

L=1.0; k=1.0; 
q0 =-1.0; fL=0.0;

ne=16; ratio=5.0;
ne=4; ratio=2.0;

%---------------
% grid generation
%---------------

xe = elm_line1(0,L,ne,ratio);

%------------------------------
% number of unique global nodes
%------------------------------

ng = 2*ne+1;

%--------------------------
% generate the global nodes
%--------------------------

Ic = 0;  % counter

for i=1:ne
 Ic = Ic+1; xg(Ic) = xe(i);
 Ic = Ic+1; xg(Ic) = 0.5*(xe(i)+xe(i+1));
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

[at,bt,ct,b] = sdqc_cnd_sys (ne,xe,q0,fL,k,s);

%--------------
% linear solver
%--------------

sol = thomas (ne,at,bt,ct,b);

%---------------------
% recover the solution
%---------------------

for i=1:ne
 f(2*i-1)=sol(i);
end
f(ng)=fL;

%-----
% plot
%-----

for i=1:ne+1
 yplot(i)=f(2*i-1);
end

plot(xe, yplot,'k-o');
for i=1:ne+1
 plot([xe(i),xe(i)], [0, 2.5],'k:');
end

%-----
% done
%-----


clear all
close all

%==============================================
% FSELIB
%
% CODE sdl_robin
%
% Code for steady one-dimensional diffusion
% with linear elements
% using the Dirichlet/Robin boundary condition
%==============================================

%-----------
% input data
%-----------

L=1.0; k=1.0; q0=-1.0; fL=0.0; finfty= 0.0;
hT=10;

ne=10; ratio=2.0;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%-------------------
% specify the source
%-------------------

for i=1:ne+1
  s(i) = 10.0*exp(-5.0*xe(i)^2/L^2);
end

%-----------------
% compact assembly
%-----------------

[at,bt,ct,b] = sdl_sys_robin (ne,xe,q0,fL,k,s,finfty,hT);

%--------------
% linear solver
%--------------

f = thomas (ne,at,bt,ct,b);

f(ne+1) = fL;

%-----
% plot
%-----

plot(xe, f,'-o');
hold on
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)

%-----
% done
%-----

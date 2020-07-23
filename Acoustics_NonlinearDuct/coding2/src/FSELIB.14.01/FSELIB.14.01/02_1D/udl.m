clear all
close all

%==============================================
% FSELIB
%
% Code udl
%
% Finite-element code for unsteady diffusion
% with linear elements
%
% SYMBOLS:
%
% ne:  number of elements
%==============================================

%-----
% prepare
%-----

figure(1)
hold on;
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)

%====================

%-----------
% input data
%-----------

L=1.0; k=1.0; rho=1.0; cp=1.0; q0=-1.0; fL=0.1;

Dt = 0.0050;
Dt = 0.00503;

Dt = 0.010;
Dt = 0.0050;

Dt = 0.0017;
Dt = 0.0016;

nsteps = 2000; nplot=100;

ne=10; ratio=1.0;

%--------
% prepare
%--------

kappa = k/(rho*cp);

%----------------
% grid generation 
%----------------

xe = elm_line1 (0,L,ne,ratio);

%------------------
% initial condition
%------------------

for i=1:ne+1
  f(1,i) = fL;
end

%-------
% source
%-------

for i=1:ne+1
  s(i) = 10.0*exp(-5.0*xe(i)^2/L^2);
end

%-----
% plot
%-----

plot(xe, f,'k+');

%==============
% time stepping
%==============

icount = 0;    % step counter

for irun=1:nsteps

%--------------------------------
% generate the tridiagonal system
%--------------------------------

  [at,bt,ct,rhs] = udl_sys (ne,xe,q0,f,fL,k,kappa,s,Dt);
% [at,bt,ct,rhs] = udl_sys_lump (ne,xe,q0,f,fL,k,kappa,s,Dt);
% [at,bt,ct,rhs] = udl_sys_cn (ne,xe,q0,f,fL,k,kappa,s,Dt);

%--------------
% linear solver
%--------------

sol = thomas (ne,at,bt,ct,rhs);

f = [sol fL];     % include the value at the right end

%---------
% plotting
%---------

if(icount==nplot)
 plot(xe, f,'-ko'); icount = 0;
end

icount = icount+1;

end

%=====================
% end of time stepping
%=====================

%------
% done
%-----

clear all
close all

%==============================================
% Code uds_lob_cn
%
% Spectral-element code for unsteady diffusion
% with the Crank-Nicolson method
%
% SYMBOLS:
%
% ne:  number of elements
% np(l):  polynomial order over the lth element
% ng:  number of distinct global nodes
% c:   connectivity matrix
%==============================================

%-----------
% input data
%-----------

k=1.0; rho=1.0; cp=1.0;
q0=-1.0; fL=0.0;
L=1.0;

ne=3; ratio=5.0;  % number of elements, stretch ratio
np(1)=4;   % element nodes
np(2)=2;
np(3)=4;  

dt=0.01;

nsteps = 150;
nplot = 5;

%--------
% prepare
%--------

kappa = k/(rho*cp);

figure(1)
hold on;
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)
box on

%------------------------
% element node generation 
%------------------------

[xe,xen,xien,xg,c,ng] = discr_lob (0,L,ne,ratio,np);

%------------------
% initial condition
%------------------

for i=1:ng
  f(1,i) = fL;
end

%-------
% source
%-------

for i=1:ng
  s(i) = 10.0*exp(-5.0*xg(i)^2/L^2);
end

%-----
% plot
%-----

plot(xg, f,'-ko');

for i=1:ne+1
 plot([xe(i),xe(i)], [0, 2.5],'k:');
end

%-----------------------------------
% global diffusion and mass matrices
%-----------------------------------

[gdm,gmm] = gdmm_lob (ne,xe,np,ng,c);

%---------------------------------
% generate the vector b on the rhs
%--------------------------------

b=zeros(1,ng-1); b(1)=q0/k;

for ip=1:ng-1
  for jp=1:ng
    b(ip) = b(ip) + gmm(ip,jp)*s(jp)/k;
  end
end

m = np(ne);
for i=1:m
  i1 = c(ne,i);
  b(i1) = b(i1) - gdm(i1,ng)*fL;
end

%-------------------------------
% generate the matrix on the LHS
%-------------------------------

lhs = zeros(ng-1,ng-1);

for i=1:ng-1
  for j=1:ng-1
    lhs(i,j) = gmm(i,j) + 0.5*dt*kappa*gdm(i,j);
  end
end

%====================
% BEGIN TIME STEPPING
%====================

icount = 0; % step counter

for irun=1:nsteps

%-----------------
% generate the RHS
%-----------------

rhs = zeros(1,ng-1);
for ip=1:ng-1
  rhs(ip) = dt*kappa*b(ip);
  for jp=1:ng-1
    rhs(ip) = rhs(ip) + (gmm(ip,jp)  ...
            - 0.5*dt*kappa*gdm(ip,jp))*f(jp);
  end
end

%--------------
% linear solver
%--------------

sol = rhs/lhs';    % solve the (ng-1)x(ng-1) linear system
f = [sol fL];     % add the value at the right end

%---------
% plotting
%---------

if(icount==nplot)
 plot(xg, f,'-ko');icount = 0;
end

icount = icount+1;
end

%=====================
% end of time stepping
%=====================

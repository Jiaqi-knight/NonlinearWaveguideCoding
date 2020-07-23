clear all
close all

%=================================
% FSELIB
%
% Code beam
%
% Cantilever beam bending
% with cubic Hermitian elements
% subject to nodal forcing
%=================================

%-----------
% input data
%-----------

E = 1.0; I = 1.0;
L = 1.0;
ne = 16; ratio = 1.0;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

ng = ne+1;  % number of global nodes

%-----------------------
% specify the nodal load
%-----------------------

for i=1:ng
 F(i) = 1.0/ne;
% F(i) = 0.0;
end
% F(ng)=0.5;

%----------------
% stifness matrix
%----------------

stiff = beam_sys(ne,xe);

stiff = E*I*stiff;

%----------------------------
% preliminary right-hand side
%----------------------------

for i=1:ng
 b(2*i-1) = - F(i);
 b(2*i) = 0.0;
end

%-------------------------------------------
% implement the left-end boundary conditions
%-------------------------------------------

for i=1:2*ng
 stiff(1,i) = 0.0;
 stiff(2,i) = 0.0;
 stiff(i,1) = 0.0;
 stiff(i,2) = 0.0;
end

stiff(1,1) = 1.0;
stiff(2,2) = 1.0;

b(1) = 0.0; b(2) = 0.0;

%--------------
% linear solver
%--------------

sol = b/stiff';

%------------------------
% extract the deflections
%------------------------

for i=1:ng
 v(i) = sol(2*i-1);
end

%----------------
% prepare to plot
%----------------

figure(1)
hold on;
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)
axis equal
axis([0 1.1*L -0.2*L 0.10*L])
box on

%---------------
% plot the nodes
%---------------

plot(xe, zeros(ng),'-ko');

%----------------------
% plot the nodal forces
%----------------------

for i=1:ng

 plotx(1) = xe(i); ploty(1) = 0.0;
 plotx(2) = xe(i); ploty(2) =-F(i);
 plot(plotx, ploty,'k');

 plttx(1) = xe(i); pltty(1) =-F(i);
 plot(plttx, pltty,'kv');

end

%---------------------
% plot the deflections
%---------------------

plot(xe,v,'.-k');

%-----
% done
%-----

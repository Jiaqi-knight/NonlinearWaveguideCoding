clear all
close all

%========================================
% Code discr_lob_dr
%
% Driver for discretizing the interval (0, L)
% into "ne" elements,
% and defining the Lobatto element
% interpolation nodes
%========================================

%-----------
% Input data
%-----------

L=1.0;
ne=2;      % number of elements
ratio=3.0;  % stretch ratio
np(1)=6; np(2) = 4;  % element expansion order

%-----------
% discretize
%-----------

[xe,xen,xien,xg,c,ng] = discr_lob (0,L,ne,ratio,np);

%-----
% plot
%-----

figure(1)
hold on;

yg= zeros(ng,1);
plot(xg,yg,'-kd','Markersize',5);   % global nodes
ye= zeros(ne+1,1);
plot(xe,ye,'ro','Markersize',10)    % end nodes
axis equal
axis([ 0 L -0.10*L 0.10*L])
xlabel('x','fontsize',15)
set(gca,'fontsize',15)
axis off
box on

%-----
% done
%-----

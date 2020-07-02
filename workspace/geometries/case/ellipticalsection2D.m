function[ellipse]=ellipticalsection2D(logfullfile,a,b,x0,y0,Neta1,Neta2,eta1min,eta1max,eta2min,eta2max,deltaeta1,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 15th, 2014
%    Last update: July 15th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);

c1(1,1) = x0;
c1(1,2) = y0-b;

c2(1,1) = x0+a;
c2(1,2) = y0;

c3(1,1) = x0;
c3(1,2) = y0+b;

c4(1,1) = x0+(-a);
c4(1,2) = y0;

k1 = (1:Neta1-2)';
k2 = (1:Neta2-2)';

theta1 = (0.5.*(2.*k1-1).*pi./(Neta1-2));
theta2 = (0.5.*(2.*k2-1).*pi./(Neta2-2));

xmin = 0;
xmax = a;
x = [c1(1,1);0.5*(xmin+xmax)+0.5*(xmax-xmin)*cos(wrev(theta1));c2(1,1)];
e1(:,1) = x0+x;
e1(:,2) = y0-b.*sqrt(1-(x./a).^2);

xmin = 0;
xmax = a;
x = [c2(1,1);wrev(0.5*(xmin+xmax)+0.5*(xmax-xmin)*cos(wrev(theta2)));c3(1,1)];
e4(:,1) = x0+x;
e4(:,2) = y0+b.*sqrt(1-(x./a).^2);

xmin = -a;
xmax = 0;
x = [c4(1,1);0.5*(xmin+xmax)+0.5*(xmax-xmin)*cos(wrev(theta1));c3(1,1)];
e3(:,1) = x0+x;
e3(:,2) = y0+b.*sqrt(1-(x./a).^2);

xmin = -a;
xmax = 0;
x = [c1(1,1);wrev(0.5*(xmin+xmax)+0.5*(xmax-xmin)*cos(wrev(theta2)));c4(1,1)];
e2(:,1) = x0+x;
e2(:,2) = y0-b.*sqrt(1-(x./a).^2);

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

ellipse = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

return
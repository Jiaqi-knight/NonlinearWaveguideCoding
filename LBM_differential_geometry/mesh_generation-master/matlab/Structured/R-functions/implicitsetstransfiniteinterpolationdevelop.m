%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 14th, 2014
%    Last update: July 14th, 2014
% 
%%

clear all
close all
clc

x0 = 1;
y0 = 2;
Lx = 50;
Ly = 20;

xC1 = 35;
yC1 = 3;
RC1 = 5;

xC2 = -15;
yC2 = 16;
RC2 = 3;


omegaxband = inline('(L^2-(x-x0).^2)./L^2','x','y','x0','L');
omegayband = inline('(L^2-(y-y0).^2)./L^2','x','y','y0','L');
omegacirc = inline('((x-xC).^2+(y-yC).^2-R^2)./R^2','x','y','xC','yC','R');

xval = ((x0-(Lx+0)):(2*Lx)/100:(x0+(Lx+0)))';
yval = ((y0-(Ly+0)):(2*Ly)/100:(y0+(Ly+0)))';

[X,Y] = meshgrid(xval,yval);

omega1 = omegaxband(X,Y,x0,Lx);
omega2 = omegayband(X,Y,y0,Ly);
omega3 = omegacirc(X,Y,xC1,yC1,RC1);
omega4 = omegacirc(X,Y,xC2,yC2,RC2);

alpha = 0;
m = 0;

omega = Rconjuction(Rconjuction(omega1,omega2,alpha,m),Rconjuction(omega3,omega4,alpha,m),alpha,m);
omegaminus = 0.5*(abs(omega)-omega);
[Fxomega,Fyomega] = gradient(omega,(2*Lx)/100,(2*Ly)/100);
[Fxomegaminus,Fyomegaminus] = gradient(omegaminus,(2*Lx)/100,(2*Ly)/100);

f1 = figure;
title('Transfinite interpolation over implicitly defined sets')
hold on
subplot(2,2,1)
mesh(X,Y,omega);
grid on
hold on
xlabel('x')
ylabel('y')
zlabel('\omega')
title('Function $\omega$','Interpreter','LaTex')
subplot(2,2,2)
mesh(X,Y,omegaminus);
grid on
hold on
xlabel('x')
ylabel('y')
zlabel('\omega')
title('Function $\frac{|\omega|-\omega}{2}$','Interpreter','LaTex')
subplot(2,2,3)
contour(omega,100);
hold on
quiver(Fxomega,Fyomega,'.k');
hold on
xlabel('x')
ylabel('y')
title('Contour plot of function $\omega$','Interpreter','LaTex')
subplot(2,2,4)
contour(omegaminus,100);
hold on
quiver(Fxomegaminus,Fyomegaminus,'.k');
hold on
xlabel('x')
ylabel('y')
title('Contour plot of function $\frac{|\omega|-\omega}{2}$','Interpreter','LaTex')

% g1 = inline('r.*sin(phi).*cos(theta)','r','phi','theta');
% g2 = inline('r.*sin(phi).*sin(theta)','r','phi','theta');
% g3 = inline('r.*cos(phi)','r','phi','theta');
% 
% rval = 1;
% 
% phival = (0:0.1:0.5*pi)';
% thetaval = (0:0.1:0.5*pi)';
% 
% [PHI,THETA] = meshgrid(phival,thetaval);
% 
% Xsphere = g1(rval,PHI,THETA);
% Ysphere = g2(rval,PHI,THETA);
% Zsphere = g3(rval,PHI,THETA);
% 
% f2 = figure;
% title('Surface of sphere')
% hold on
% subplot(2,3,1)
% mesh(PHI,THETA,Xsphere);
% grid on
% hold on
% xlabel('\phi')
% ylabel('\theta')
% zlabel('X')
% title('X as a function of $\phi$ and $\theta$','Interpreter','LaTex')
% subplot(2,3,2)
% mesh(PHI,THETA,Ysphere);
% grid on
% hold on
% xlabel('\phi')
% ylabel('\theta')
% zlabel('Y')
% title('Y as a function of $\phi$ and $\theta$','Interpreter','LaTex')
% subplot(2,3,3)
% mesh(PHI,THETA,Zsphere);
% grid on
% hold on
% xlabel('\phi')
% ylabel('\theta')
% zlabel('Z')
% title('Z as a function of $\phi$ and $\theta$','Interpreter','LaTex')
% subplot(2,3,4)
% contour(Xsphere,100);
% hold on
% xlabel('\phi')
% ylabel('\theta')
% title('Contour plot of X as a function of $\phi$ and $\theta$','Interpreter','LaTex')
% subplot(2,3,5)
% contour(Ysphere,100);
% hold on
% xlabel('\phi')
% ylabel('\theta')
% title('Contour plot of Y as a function of $\phi$ and $\theta$','Interpreter','LaTex')
% subplot(2,3,6)
% contour(Zsphere,100);
% hold on
% xlabel('\phi')
% ylabel('\theta')
% title('Contour plot of Z as a function of $\phi$ and $\theta$','Interpreter','LaTex')
%%
%==============================================================================
% Copyright (c) 2016 Universit? de Lorraine & Lule? tekniska universitet
% Author: Luca Di Stasio <luca.distasio@gmail.com>
%                        <luca.distasio@ingpec.eu>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 
% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% Neither the name of the Universit? de Lorraine or Lule? tekniska universitet
% nor the names of its contributors may be used to endorse or promote products
% derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================
%
%  DESCRIPTION
%  
%  A function to perform 
%
%  Input: 
%  Output: 
%
%%

clear all
close all
clc
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
% subfunction_path4=genpath('C:\Users\wjq\Desktop\workspace\geometry-master\geometry')
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
% addpath(subfunction_path4);

formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta1min = 0;
Neta1 = 50;
deltaeta1 = 1;
eta1max = eta1min + Neta1*deltaeta1;
eta2min = 0;
Neta2 = 60;
deltaeta2 = 1;
eta2max = eta2min + Neta2*deltaeta2;

N = Neta1*Neta2;

xmin = 0;
xmax = 5;
Nx = Neta1;
ymin = -2.5;
ymax = 3;
Ny = Neta2;

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

%                 e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
%                 e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
%                 e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
%                 e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
%                 c1 = c1(xi_min,eta_min) = (x1,y1)
%                 c2 = c2(xi_max,eta_min) = (x2,y2)
%                 c3 = c3(xi_max,eta_max) = (x3,y3)
%                 c4 = c4(xi_min,eta_max) = (x4,y4)
%                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

c1(1,1:2) = [1 1];
c2(1,1:2) = [10 1];
c3(1,1:2) = [10 5];
c4(1,1:2) = [1 5];

xcm = 0.5*(c3(1,1)+c4(1,1));
ycm = c3(1,2);
r = xcm - c4(1,1);
theta = (pi:(0-pi)/(Neta1-1):0)';
e3(:,1) = xcm*ones(Neta1,1)+r.*(1-0.5.*abs(sin(theta))).*cos(theta);
e3(:,2) = ycm*ones(Neta1,1)+r.*(1-0.5.*abs(sin(theta))).*sin(theta);

xcm = 0.5*(c1(1,1)+c2(1,1));
ycm = c1(1,2);
r = xcm - c1(1,1);
theta = (pi:(0-pi)/(Neta1-1):0)';
e1(:,1) = xcm*ones(Neta1,1)+r.*(1-0.5.*abs(sin(theta))).*cos(theta);
e1(:,2) = ycm*ones(Neta1,1)-r.*(1-0.5.*abs(sin(theta))).*sin(theta);

ycm = 0.5*(c1(1,2)+c4(1,2));
xcm = c1(1,1);
r = ycm - c1(1,2);
theta = (pi:(0-pi)/(Neta2-1):0)';
e2(:,1) = xcm*ones(Neta2,1)+r.*(1-0.5.*abs(sin(theta))).*sin(theta);
e2(:,2) = ycm*ones(Neta2,1)+r.*(1-0.5.*abs(sin(theta))).*cos(theta);

ycm = 0.5*(c2(1,2)+c3(1,2));
xcm = c2(1,1);
r = ycm - c2(1,2);
theta = (pi:(0-pi)/(Neta2-1):0)';
e4(:,1) = xcm*ones(Neta2,1)-r.*(1-0.5.*abs(sin(theta))).*sin(theta);
e4(:,2) = ycm*ones(Neta2,1)+r.*(1-0.5.*abs(sin(theta))).*cos(theta);


f = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
plot(e1(:,1),e1(:,2),'r*')
hold on
plot(e2(:,1),e2(:,2),'r*')
hold on
plot(e3(:,1),e3(:,2),'r*')
hold on
plot(e4(:,1),e4(:,2),'r*')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
axis([0 11 -5 11])

lattice = transfiniteinterpolation2D(logfullfile,N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice(:,1),lattice(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice(:,3),lattice(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
axis(fphys,[0 11 -5 11])

clearvars -except logfullfile % 只保留变量x , y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta1min = 0;
Neta1 = 50;
deltaeta1 = 1;
eta1max = eta1min + Neta1*deltaeta1;
eta2min = 0;
Neta2 = 60;
deltaeta2 = 1;
eta2max = eta2min + Neta2*deltaeta2;

N = Neta1*Neta2;

xmin = 0;
xmax = 5;
Nx = Neta1;
ymin = -2.5;
ymax = 3;
Ny = Neta2;

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

%                 e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
%                 e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
%                 e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
%                 e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
%                 c1 = c1(xi_min,eta_min) = (x1,y1)
%                 c2 = c2(xi_max,eta_min) = (x2,y2)
%                 c3 = c3(xi_max,eta_max) = (x3,y3)
%                 c4 = c4(xi_min,eta_max) = (x4,y4)
%                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

center = [5 5];
r = 10;

x0 = center(1);
y0 = center(2);

l = 5;

c1(1,1:2) = [x0-l y0+l];
c2(1,1:2) = [x0+l y0+l];
c3(1,1:2) = [x0+r/sqrt(2) y0+r/sqrt(2)];
c4(1,1:2) = [x0-r/sqrt(2) y0+r/sqrt(2)];

xs = (c1(1):2*l/(Neta1-1):c2(1))';
e1(:,1) = xs;
e1(:,2) = c1(2)*ones(length(xs),1);

e4(:,1) = (c2(1):(c3(1)-c2(1))/(Neta2-1):c3(1))';
e4(:,2) = (c2(2):(c3(2)-c2(2))/(Neta2-1):c3(2))';

theta = (0.75*pi:-0.5*pi/(Neta1-1):0.25*pi)';
e3(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
e3(:,2) = y0*ones(Neta1,1)+r.*sin(theta);

%xs = (c4(1):(c3(1)-c4(1))/(Neta1-1):c3(1))';
%e3(:,1) = xs;
%e3(:,2) = c3(2)*ones(length(xs),1);

e2(:,1) = (c1(1):(c4(1)-c1(1))/(Neta2-1):c4(1))';
e2(:,2) = (c1(2):(c4(2)-c1(2))/(Neta2-1):c4(2))';


f = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
plot(e1(:,1),e1(:,2),'r*')
hold on
plot(e2(:,1),e2(:,2),'r*')
hold on
plot(e3(:,1),e3(:,2),'r*')
hold on
plot(e4(:,1),e4(:,2),'r*')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
%axis([0 11 -5 11])

lattice = transfiniteinterpolation2D(logfullfile,N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice(:,1),lattice(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice(:,3),lattice(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta1min = 0;
Neta1 = 50;
deltaeta1 = 1;
eta1max = eta1min + Neta1*deltaeta1;
eta2min = 0;
Neta2 = 60;
deltaeta2 = 1;
eta2max = eta2min + Neta2*deltaeta2;

N = Neta1*Neta2;

xmin = 0;
xmax = 5;
Nx = Neta1;
ymin = -2.5;
ymax = 3;
Ny = Neta2;

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

% %                 e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
% %                 e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
% %                 e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
% %                 e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
% %                 c1 = c1(xi_min,eta_min) = (x1,y1)
% %                 c2 = c2(xi_max,eta_min) = (x2,y2)
% %                 c3 = c3(xi_max,eta_max) = (x3,y3)
% %                 c4 = c4(xi_min,eta_max) = (x4,y4)
% %                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
% %                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
% %                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
% %                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
% %                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
% %                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
% %                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
% %                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

center = [5 5];
r = 10;

x0 = center(1);
y0 = center(2);

l = 5;

c1(1,1:2) = [x0-r/sqrt(2) y0-r/sqrt(2)];
c2(1,1:2) = [x0+r/sqrt(2) y0-r/sqrt(2)];
c3(1,1:2) = [x0+r/sqrt(2) y0+r/sqrt(2)];
c4(1,1:2) = [x0-r/sqrt(2) y0+r/sqrt(2)];

theta = (5*pi/4:0.5*pi/(Neta1-1):7*pi/4)';
e1(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
e1(:,2) = y0*ones(Neta1,1)+r.*sin(theta);

theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e2(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
e2(:,2) = y0*ones(Neta2,1)+r.*sin(theta);

theta = (3*pi/4:-0.5*pi/(Neta1-1):0.25*pi)';
e3(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
e3(:,2) = y0*ones(Neta1,1)+r.*sin(theta);

theta = (-0.25*pi:0.5*pi/(Neta2-1):0.25*pi)';
e4(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
e4(:,2) = y0*ones(Neta2,1)+r.*sin(theta);


f = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
plot(e1(:,1),e1(:,2),'r*')
hold on
plot(e2(:,1),e2(:,2),'r*')
hold on
plot(e3(:,1),e3(:,2),'r*')
hold on
plot(e4(:,1),e4(:,2),'r*')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
%axis([0 11 -5 11])

lattice = transfiniteinterpolation2D(N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice(:,1),lattice(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice(:,3),lattice(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

figure();
plot(lattice(indicesC1,3),lattice(indicesC1,4),'.')
hold on
plot(lattice(indicesC2,3),lattice(indicesC2,4),'.')
hold on
plot(lattice(indicesC3,3),lattice(indicesC3,4),'.')
hold on
plot(lattice(indicesC4,3),lattice(indicesC4,4),'.')
hold on
plot(lattice(indicesE1,3),lattice(indicesE1,4),'.')
hold on
plot(lattice(indicesE2,3),lattice(indicesE2,4),'.')
hold on
plot(lattice(indicesE3,3),lattice(indicesE3,4),'.')
hold on
plot(lattice(indicesE4,3),lattice(indicesE4,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - transfinite interpolation')
axis equal

flagperiodicity = 0;
periodicity = 0;
flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;
deltaq = [deltaeta1 deltaeta2];
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(N,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

itmax = 10;
tol = 10^-4;

lattice10 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,1); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice10(:,1),lattice10(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice10(:,5),lattice10(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 10 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])


itmax = 50;
tol = 10^-4;

lattice50 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice50(:,1),lattice50(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice50(:,5),lattice50(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 50 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])


itmax = 100;
tol = 10^-4;

lattice100 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice100(:,1),lattice100(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice100(:,5),lattice100(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 100 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 1000;
tol = 10^-4;

lattice1000 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice1000(:,1),lattice1000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice1000(:,5),lattice1000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 5000;
tol = 10^-4;

lattice5000 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice5000(:,1),lattice5000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice5000(:,5),lattice5000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 5000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])


itmax = 10000;
tol = 10^-4;

lattice10000 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice10000(:,1),lattice10000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice10000(:,5),lattice10000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 10000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 50000;
tol = 10^-4;

lattice50000 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice50000(:,1),lattice50000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice50000(:,5),lattice50000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 50000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 100000;
tol = 10^-4;

lattice100000 = sparseellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

scrsz = get(0,'ScreenSize');

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice100000(:,1),lattice100000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice100000(:,5),lattice100000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 100000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% eta1min = 0;
% Neta1 = 50;
% deltaeta1 = 1;
% eta1max = eta1min + Neta1*deltaeta1;
% eta2min = 0;
% Neta2 = 60;
% deltaeta2 = 1;
% eta2max = eta2min + Neta2*deltaeta2;
% 
% N = Neta1*Neta2;
% 
% xmin = 0;
% xmax = 5;
% Nx = Neta1;
% ymin = -2.5;
% ymax = 3;
% Ny = Neta2;
% 
% lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
% 
% [indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);
% 
% e1 = zeros(Neta1,6);
% e2 = zeros(Neta2,6);
% e3 = zeros(Neta1,6);
% e4 = zeros(Neta2,6);
% c1 = zeros(1,8);
% c2 = zeros(1,8);
% c3 = zeros(1,8);
% c4 = zeros(1,8);
% 
% center = [5 5];
% r = 10;
% 
% x0 = center(1);
% y0 = center(2);
% 
% l = 5;
% 
% c1(1,1:2) = [x0-r/sqrt(2) y0-r/sqrt(2)];
% eta = 5*pi/4;
% xi  = r;
% c1(1,3:8) = [cos(eta) -xi*sin(eta) sin(eta) xi*cos(eta) -sin(eta) cos(eta)];
% c2(1,1:2) = [x0+r/sqrt(2) y0-r/sqrt(2)];
% eta = 7*pi/4;
% xi  = r;
% c2(1,3:8) = [cos(eta) -xi*sin(eta) sin(eta) xi*cos(eta) -sin(eta) cos(eta)];
% c3(1,1:2) = [x0+r/sqrt(2) y0+r/sqrt(2)];
% eta = pi/4;
% xi  = r;
% c3(1,3:8) = [cos(eta) -xi*sin(eta) sin(eta) xi*cos(eta) -sin(eta) cos(eta)];
% c4(1,1:2) = [x0-r/sqrt(2) y0+r/sqrt(2)];
% eta = 3*pi/4;
% xi  = r;
% c4(1,3:8) = [cos(eta) -xi*sin(eta) sin(eta) xi*cos(eta) -sin(eta) cos(eta)];
% 
% theta = (5*pi/4:0.5*pi/(Neta1-1):7*pi/4)';
% e1(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
% e1(:,2) = y0*ones(Neta1,1)+r.*sin(theta);
% xi  = r;
% e1(:,3:6) = [cos(theta) -xi*sin(theta) sin(theta) xi*cos(theta)];
% 
% theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
% e2(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
% e2(:,2) = y0*ones(Neta2,1)+r.*sin(theta);
% xi  = r;
% e2(:,3:6) = [cos(theta) -xi*sin(theta) sin(theta) xi*cos(theta)];
% 
% theta = (3*pi/4:-0.5*pi/(Neta1-1):0.25*pi)';
% e3(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
% e3(:,2) = y0*ones(Neta1,1)+r.*sin(theta);
% xi  = r;
% e3(:,3:6) = [cos(theta) -xi*sin(theta) sin(theta) xi*cos(theta)];
% 
% theta = (-0.25*pi:0.5*pi/(Neta2-1):0.25*pi)';
% e4(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
% e4(:,2) = y0*ones(Neta2,1)+r.*sin(theta);
% xi  = r;
% e4(:,3:6) = [cos(theta) -xi*sin(theta) sin(theta) xi*cos(theta)];
% 
% f = figure();
% plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
% hold on
% plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
% hold on
% plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
% hold on
% plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
% hold on
% plot(e1(:,1),e1(:,2),'r*')
% hold on
% plot(e2(:,1),e2(:,2),'r*')
% hold on
% plot(e3(:,1),e3(:,2),'r*')
% hold on
% plot(e4(:,1),e4(:,2),'r*')
% hold on
% grid on
% xlabel('x')
% ylabel('y')
% title('Mesh in fluid domain (physical space)')
% axis equal
% %axis([0 11 -5 11])
% 
% lattice = transfiniteinterpolation2D(N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,2,e1,e2,e3,e4,c1,c2,c3,c4);
% 
% scrsz = get(0,'ScreenSize');
% 
% f1 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
% title('Computational and physical domains')
% hold on
% 
% subplot(1,2,1);
% subplot(1,2,2);
% 
% hsubfigs1 = get(f1,'Children');
% 
% fcomp1 = hsubfigs1(2);
% fphys1 = hsubfigs1(1);
% 
% subplot(fcomp1);
% plot(lattice(:,1),lattice(:,2),'.')
% hold on
% grid on
% xlabel('$\xi_{1}$','Interpreter','LaTex')
% ylabel('$\xi_{2}$','Interpreter','LaTex')
% title('Mesh in lattice domain (computational space)')
% 
% subplot(fphys1);
% plot(lattice(:,5),lattice(:,6),'.')
% hold on
% grid on
% xlabel('x')
% ylabel('y')
% title('Mesh in fluid domain (physical space)')
% 
% axis([fcomp fphys],'equal')
% axis([fcomp],[0 100 0 100])
% % axis([fphys],[0 11 -5 11])
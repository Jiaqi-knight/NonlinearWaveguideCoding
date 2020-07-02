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

N = 100;

a = 10;
b = 3;
x0 = 0;
y0 = 0;

ellipse = ellipse2D(a,b,0,x0,y0,N);

f0 = figure;
plot(ellipse(:,1),ellipse(:,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Ellipse')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 10;
b = 3;
x0 = 0;
y0 = 0;

ellipse = ellipse2D(a,b,0,x0,y0,N);

f1 = figure;
plot(ellipse(:,1),ellipse(:,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Ellipse')
axis equal

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

%ellipse = ellipticalsection2D(a,b,x0,y0,Neta1,Neta2,eta1min,eta1max,eta2min,eta2max,deltaeta1,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);

c = sqrt(a^2-b^2);

x = -c;
c1(1,1) = x0-c;
c1(1,2) = y0-b.*sqrt(1-(x./a).^2);

x = c;
c2(1,1) = x0+c;
c2(1,2) = y0-b.*sqrt(1-(x./a).^2);

x = c;
c3(1,1) = x0+c;
c3(1,2) = y0+b.*sqrt(1-(x./a).^2);

x = -c;
c4(1,1) = x0-c;
c4(1,2) = y0+b.*sqrt(1-(x./a).^2);

x = (-c:2*c/(Neta1-1):c)';
e1(:,1) = x0+x;
e1(:,2) = y0-b.*sqrt(1-(x./a).^2);

delta = (a-c)/(Neta2/2-1);
x1 = (c:(a-delta-c)/(Neta2/2-1):a-delta)';
x2 = (a-delta:-(a-delta-c)/(Neta2/2-1):c)';
x = [x1;x2];
e4(:,1) = x0+x;
e4(:,2) = y0+[-b.*sqrt(1-(x1./a).^2);b.*sqrt(1-(x2./a).^2)];

x = (-c:2*c/(Neta1-1):c)';
e3(:,1) = x0+x;
e3(:,2) = y0+b.*sqrt(1-(x./a).^2);

delta = (a-c)/(Neta2/2-1);
x1 = (-c:-(a-delta-c)/(Neta2/2-1):-a+delta)';
x2 = (-a+delta:(a-delta-c)/(Neta2/2-1):-c)';
x = [x1;x2];
e2(:,1) = x0+x;
e2(:,2) = y0+[-b.*sqrt(1-(x1./a).^2);b.*sqrt(1-(x2./a).^2)];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

ellipse = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

f2 = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
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

f3 = figure;
plot(ellipse(:,3),ellipse(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 10;
b = 3;
x0 = 0;
y0 = 0;

ellipse = ellipse2D(a,b,0,x0,y0,N);

f1 = figure;
plot(ellipse(:,1),ellipse(:,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Ellipse')
axis equal

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

%ellipse = ellipticalsection2D(a,b,x0,y0,Neta1,Neta2,eta1min,eta1max,eta2min,eta2max,deltaeta1,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);

c = sqrt(a^2-b^2);

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

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,Nx,Ny);

ellipse = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

f2 = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
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

f3 = figure;
plot(ellipse(:,3),ellipse(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

figure();
plot(ellipse(indicesC1,3),ellipse(indicesC1,4),'.')
hold on
plot(ellipse(indicesC2,3),ellipse(indicesC2,4),'.')
hold on
plot(ellipse(indicesC3,3),ellipse(indicesC3,4),'.')
hold on
plot(ellipse(indicesC4,3),ellipse(indicesC4,4),'.')
hold on
plot(ellipse(indicesE1,3),ellipse(indicesE1,4),'.')
hold on
plot(ellipse(indicesE2,3),ellipse(indicesE2,4),'.')
hold on
plot(ellipse(indicesE3,3),ellipse(indicesE3,4),'.')
hold on
plot(ellipse(indicesE4,3),ellipse(indicesE4,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - transfinite interpolation')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 10;
b = 3;
x0 = 0;
y0 = 0;

ellipse = ellipse2D(a,b,0,x0,y0,N);

f1 = figure;
plot(ellipse(:,1),ellipse(:,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Ellipse')
axis equal

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

ellipse = ellipticalsection2D(logfullfile,a,b,x0,y0,Neta1,Neta2,eta1min,eta1max,eta2min,eta2max,deltaeta1,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

f2 = figure();
plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
hold on
plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
hold on
plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
hold on
plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal
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

f3 = figure;
plot(ellipse(:,3),ellipse(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];

Neta1 = 30;
eta1min = 0;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

Neta2 = 25;
eta2min = 0;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

Neta3 = 35;
eta3min = 0;
deltaeta3 = 1;
eta3max = eta3min + (Neta3-1)*deltaeta3;

Nx = Neta1;
xmin = eta1min;
xmax = eta1max;

Ny = Neta2;
ymin = eta2min;
ymax = eta2max;

Nz = Neta3;
zmin = eta3min;
zmax = eta3max;

c1 = zeros(1,3);
c2 = zeros(1,3);
c3 = zeros(1,3);
c4 = zeros(1,3);
c5 = zeros(1,3);
c6 = zeros(1,3);
c7 = zeros(1,3);
c8 = zeros(1,3);

e1  = zeros(Neta1,3);
e2  = zeros(Neta2,3);
e3  = zeros(Neta1,3);
e4  = zeros(Neta2,3);
e5  = zeros(Neta3,3);
e6  = zeros(Neta3,3);
e7  = zeros(Neta3,3);
e8  = zeros(Neta3,3);
e9  = zeros(Neta1,3);
e10 = zeros(Neta2,3);
e11 = zeros(Neta1,3);
e12 = zeros(Neta2,3);

f1 = zeros(Neta1*Neta2,3);
f2 = zeros(Neta1*Neta3,3);
f3 = zeros(Neta2*Neta3,3);
f4 = zeros(Neta1*Neta3,3);
f5 = zeros(Neta2*Neta3,3);
f6 = zeros(Neta1*Neta2,3);

R = 20;

xfunc = @(r,phi,theta) r.*sin(phi).*cos(theta);
yfunc = @(r,phi,theta) r.*sin(phi).*sin(theta);
zfunc = @(r,phi,theta) r.*cos(phi);

c1(1,:) = [xfunc(R,0.75*pi,5*pi/4) yfunc(R,0.75*pi,5*pi/4) zfunc(R,0.75*pi,5*pi/4)];
c2(1,:) = [xfunc(R,0.75*pi,7*pi/4) yfunc(R,0.75*pi,7*pi/4) zfunc(R,0.75*pi,7*pi/4)];
c3(1,:) = [xfunc(R,0.75*pi,pi/4)   yfunc(R,0.75*pi,pi/4)   zfunc(R,0.75*pi,pi/4)];
c4(1,:) = [xfunc(R,0.75*pi,3*pi/4) yfunc(R,0.75*pi,3*pi/4) zfunc(R,0.75*pi,3*pi/4)];
c5(1,:) = [xfunc(R,0.25*pi,5*pi/4) yfunc(R,0.25*pi,5*pi/4) zfunc(R,0.25*pi,5*pi/4)];
c6(1,:) = [xfunc(R,0.25*pi,7*pi/4) yfunc(R,0.25*pi,7*pi/4) zfunc(R,0.25*pi,7*pi/4)];
c7(1,:) = [xfunc(R,0.25*pi,pi/4)   yfunc(R,0.25*pi,pi/4)   zfunc(R,0.25*pi,pi/4)];
c8(1,:) = [xfunc(R,0.25*pi,3*pi/4) yfunc(R,0.25*pi,3*pi/4) zfunc(R,0.25*pi,3*pi/4)];

xs = (c1(1,1):(c2(1,1)-c1(1,1))/(Neta1-1):c2(1,1))';
e1(:,1:3) = [xs c1(1,2)*ones(length(xs),1) c1(1,3)*ones(length(xs),1)];

phi = 0.75*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e2(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% e2(:,1:3) = [c2(1,1)*ones(length(ys),1) ys c2(1,3)*ones(length(ys),1)];

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e3(:,1:3) = [xs c4(1,2)*ones(length(xs),1) c4(1,3)*ones(length(xs),1)];

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
e4(:,1:3) = [c1(1,1)*ones(length(ys),1) ys c1(1,3)*ones(length(ys),1)];

zs = (c1(1,3):(c5(1,3)-c1(1,3))/(Neta3-1):c5(1,3))';
e5(:,1:3) = [c1(1,1)*ones(length(zs),1) c1(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e6(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c2(1,1) c2(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c2(1,3):(c6(1,3)-c2(1,3))/(Neta3-1):c6(1,3))';
% e6(:,1:3) = [c2(1,1)*ones(length(zs),1) c2(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e7(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c3(1,1) c3(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c3(1,3):(c7(1,3)-c3(1,3))/(Neta3-1):c7(1,3))';
% e7(:,1:3) = [c3(1,1)*ones(length(zs),1) c3(1,2)*ones(length(zs),1) zs];

zs = (c4(1,3):(c8(1,3)-c4(1,3))/(Neta3-1):c8(1,3))';
e8(:,1:3) = [c4(1,1)*ones(length(zs),1) c4(1,2)*ones(length(zs),1) zs];

xs = (c5(1,1):(c6(1,1)-c5(1,1))/(Neta1-1):c6(1,1))';
e9(:,1:3) = [xs c5(1,2)*ones(length(xs),1) c5(1,3)*ones(length(xs),1)];

phi = 0.25*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e10(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c6(1,2):(c7(1,2)-c6(1,2))/(Neta2-1):c7(1,2))';
% e10(:,1:3) = [c6(1,1)*ones(length(ys),1) ys c6(1,3)*ones(length(ys),1)];

xs = (c8(1,1):(c7(1,1)-c8(1,1))/(Neta1-1):c7(1,1))';
e11(:,1:3) = [xs c8(1,2)*ones(length(xs),1) c8(1,3)*ones(length(xs),1)];

ys = (c5(1,2):(c8(1,2)-c5(1,2))/(Neta2-1):c8(1,2))';
e12(:,1:3) = [c5(1,1)*ones(length(ys),1) ys c5(1,3)*ones(length(ys),1)];

p1 = 0.999;
p2 = 0.001;


lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e4,e3,e2,c1,c2,c3,c4);
f1(:,1:3) = [lattice(:,5) lattice(:,6) c1(1,3)*ones(length(lattice),1)];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e1(:,1) e1(:,3)],[e5(:,1) e5(:,3)],[e9(:,1) e9(:,3)],[e6(:,1) e6(:,3)],[c1(1,1) c1(1,3)],[c2(1,1) c2(1,3)],[c6(1,1) c6(1,3)],[c5(1,1) c5(1,3)]);
f2(:,1:3) = [lattice(:,5) c1(1,2)*ones(length(lattice),1) lattice(:,6)];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
for i=1:Neta3
    f3((i-1)*Neta2+1:i*Neta2,1:3) = [e2(:,1)+e6(i,1)-c2(1,1) e2(:,2) zfunc(R,phi(i),0)*ones(length(e2),1)];
end

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% deltaz = (c6(1,3)-c2(1,3))/(Neta3-1);
% for i=1:Neta3
%     f3((i-1)*Neta2+1:i*Neta2,1:3) = [c2(1,1)*ones(length(ys),1) ys (c2(1,3)+deltaz*(i-1))*ones(length(ys),1)];
% end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e3(:,1) e3(:,3)],[e8(:,1) e8(:,3)],[e11(:,1) e11(:,3)],[e7(:,1) e7(:,3)],[c4(1,1) c4(1,3)],[c3(1,1) c3(1,3)],[c7(1,1) c7(1,3)],[c8(1,1) c8(1,3)]);
f4(:,1:3) = [lattice(:,5) c3(1,2)*ones(length(lattice),1) lattice(:,6)];

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
deltaz = (c5(1,3)-c1(1,3))/(Neta3-1);
for i=1:Neta3
    f5((i-1)*Neta2+1:i*Neta2,1:3) = [c1(1,1)*ones(length(ys),1) ys (c1(1,3)+deltaz*(i-1))*ones(length(ys),1)];
end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e9,e12,e11,e10,c5,c6,c7,c8);
f6(:,1:3) = [lattice(:,5) lattice(:,6) c5(1,3)*ones(length(lattice),1)];

fig0 = figure();
plot3(c1(1,1),c1(1,2),c1(1,3),'kx')
hold on
plot3(c2(1,1),c2(1,2),c2(1,3),'kx')
hold on
plot3(c3(1,1),c3(1,2),c3(1,3),'kx')
hold on
plot3(c4(1,1),c4(1,2),c4(1,3),'kx')
hold on
plot3(c5(1,1),c5(1,2),c5(1,3),'kx')
hold on
plot3(c6(1,1),c6(1,2),c6(1,3),'kx')
hold on
plot3(c7(1,1),c7(1,2),c7(1,3),'kx')
hold on
plot3(c8(1,1),c8(1,2),c8(1,3),'kx')
hold on
plot3(e1(:,1),e1(:,2),e1(:,3),'r*')
hold on
plot3(e2(:,1),e2(:,2),e2(:,3),'r*')
hold on
plot3(e3(:,1),e3(:,2),e3(:,3),'r*')
hold on
plot3(e4(:,1),e4(:,2),e4(:,3),'r*')
hold on
plot3(e5(:,1),e5(:,2),e5(:,3),'r*')
hold on
plot3(e6(:,1),e6(:,2),e6(:,3),'r*')
hold on
plot3(e7(:,1),e7(:,2),e7(:,3),'r*')
hold on
plot3(e8(:,1),e8(:,2),e8(:,3),'r*')
hold on
plot3(e9(:,1),e9(:,2),e9(:,3),'r*')
hold on
plot3(e10(:,1),e10(:,2),e10(:,3),'r*')
hold on
plot3(e11(:,1),e11(:,2),e11(:,3),'r*')
hold on
plot3(e12(:,1),e12(:,2),e12(:,3),'r*')
hold on
plot3(f1(:,1),f1(:,2),f1(:,3),'b.')
hold on
plot3(f2(:,1),f2(:,2),f2(:,3),'b.')
hold on
plot3(f3(:,1),f3(:,2),f3(:,3),'b.')
hold on
plot3(f4(:,1),f4(:,2),f4(:,3),'b.')
hold on
plot3(f5(:,1),f5(:,2),f5(:,3),'b.')
hold on
plot3(f6(:,1),f6(:,2),f6(:,3),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

geom = generatelattice3D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,ymin,ymax,Ny,zmin,zmax,Nz);

geom = transfiniteinterpolation3D(Neta1,Neta2,Neta3,geom(:,1:3),eta1min,eta1max,eta2min,eta2max,eta3min,eta3max,1,f1,f2,f3,f4,f5,f6,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,c1,c2,c3,c4,c5,c6,c7,c8);

fig1 = figure();
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

fig2 = figure();
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
axis equal

fig3 = figure();
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(fig3,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except logfullfile % 只保留变量x , y

Neta1 = 30;
eta1min = 0;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

Neta2 = 25;
eta2min = 0;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

Neta3 = 35;
eta3min = 0;
deltaeta3 = 1;
eta3max = eta3min + (Neta3-1)*deltaeta3;

Nx = Neta1;
xmin = eta1min;
xmax = eta1max;

Ny = Neta2;
ymin = eta2min;
ymax = eta2max;

Nz = Neta3;
zmin = eta3min;
zmax = eta3max;

c1 = zeros(1,3);
c2 = zeros(1,3);
c3 = zeros(1,3);
c4 = zeros(1,3);
c5 = zeros(1,3);
c6 = zeros(1,3);
c7 = zeros(1,3);
c8 = zeros(1,3);

e1  = zeros(Neta1,3);
e2  = zeros(Neta2,3);
e3  = zeros(Neta1,3);
e4  = zeros(Neta2,3);
e5  = zeros(Neta3,3);
e6  = zeros(Neta3,3);
e7  = zeros(Neta3,3);
e8  = zeros(Neta3,3);
e9  = zeros(Neta1,3);
e10 = zeros(Neta2,3);
e11 = zeros(Neta1,3);
e12 = zeros(Neta2,3);

f1 = zeros(Neta1*Neta2,3);
f2 = zeros(Neta1*Neta3,3);
f3 = zeros(Neta2*Neta3,3);
f4 = zeros(Neta1*Neta3,3);
f5 = zeros(Neta2*Neta3,3);
f6 = zeros(Neta1*Neta2,3);

R = 20;

xfunc = @(r,phi,theta) r.*sin(phi).*cos(theta);
yfunc = @(r,phi,theta) r.*sin(phi).*sin(theta);
zfunc = @(r,phi,theta) r.*cos(phi);

c1(1,:) = [xfunc(R,0.75*pi,5*pi/4) yfunc(R,0.75*pi,5*pi/4) zfunc(R,0.75*pi,5*pi/4)];
c2(1,:) = [xfunc(R,0.75*pi,7*pi/4) yfunc(R,0.75*pi,7*pi/4) zfunc(R,0.75*pi,7*pi/4)];
c3(1,:) = [xfunc(R,0.75*pi,pi/4)   yfunc(R,0.75*pi,pi/4)   zfunc(R,0.75*pi,pi/4)];
c4(1,:) = [xfunc(R,0.75*pi,3*pi/4) yfunc(R,0.75*pi,3*pi/4) zfunc(R,0.75*pi,3*pi/4)];
c5(1,:) = [xfunc(R,0.25*pi,5*pi/4) yfunc(R,0.25*pi,5*pi/4) zfunc(R,0.25*pi,5*pi/4)];
c6(1,:) = [xfunc(R,0.25*pi,7*pi/4) yfunc(R,0.25*pi,7*pi/4) zfunc(R,0.25*pi,7*pi/4)];
c7(1,:) = [xfunc(R,0.25*pi,pi/4)   yfunc(R,0.25*pi,pi/4)   zfunc(R,0.25*pi,pi/4)];
c8(1,:) = [xfunc(R,0.25*pi,3*pi/4) yfunc(R,0.25*pi,3*pi/4) zfunc(R,0.25*pi,3*pi/4)];

xs = (c1(1,1):(c2(1,1)-c1(1,1))/(Neta1-1):c2(1,1))';
e1(:,1:3) = [xs c1(1,2)*ones(length(xs),1) c1(1,3)*ones(length(xs),1)];

phi = 0.75*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e2(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% e2(:,1:3) = [c2(1,1)*ones(length(ys),1) ys c2(1,3)*ones(length(ys),1)];

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e3(:,1:3) = [xs c4(1,2)*ones(length(xs),1) c4(1,3)*ones(length(xs),1)];

phi = 0.75*pi;
theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e4(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
% e4(:,1:3) = [c1(1,1)*ones(length(ys),1) ys c1(1,3)*ones(length(ys),1)];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = pi;
e5(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c1(1,1) c1(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c1(1,3):(c5(1,3)-c1(1,3))/(Neta3-1):c5(1,3))';
% e5(:,1:3) = [c1(1,1)*ones(length(zs),1) c1(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e6(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c2(1,1) c2(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c2(1,3):(c6(1,3)-c2(1,3))/(Neta3-1):c6(1,3))';
% e6(:,1:3) = [c2(1,1)*ones(length(zs),1) c2(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e7(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c3(1,1) c3(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c3(1,3):(c7(1,3)-c3(1,3))/(Neta3-1):c7(1,3))';
% e7(:,1:3) = [c3(1,1)*ones(length(zs),1) c3(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = pi;
e8(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c4(1,1) c4(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c4(1,3):(c8(1,3)-c4(1,3))/(Neta3-1):c8(1,3))';
% e8(:,1:3) = [c4(1,1)*ones(length(zs),1) c4(1,2)*ones(length(zs),1) zs];

xs = (c5(1,1):(c6(1,1)-c5(1,1))/(Neta1-1):c6(1,1))';
e9(:,1:3) = [xs c5(1,2)*ones(length(xs),1) c5(1,3)*ones(length(xs),1)];

phi = 0.25*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e10(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c6(1,2):(c7(1,2)-c6(1,2))/(Neta2-1):c7(1,2))';
% e10(:,1:3) = [c6(1,1)*ones(length(ys),1) ys c6(1,3)*ones(length(ys),1)];

xs = (c8(1,1):(c7(1,1)-c8(1,1))/(Neta1-1):c7(1,1))';
e11(:,1:3) = [xs c8(1,2)*ones(length(xs),1) c8(1,3)*ones(length(xs),1)];

phi = 0.25*pi;
theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e12(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c5(1,2):(c8(1,2)-c5(1,2))/(Neta2-1):c8(1,2))';
% e12(:,1:3) = [c5(1,1)*ones(length(ys),1) ys c5(1,3)*ones(length(ys),1)];

p1 = 0.999;
p2 = 0.001;


lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e4,e3,e2,c1,c2,c3,c4);
f1(:,1:3) = [lattice(:,5) lattice(:,6) c1(1,3)*ones(length(lattice),1)];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e1(:,1) e1(:,3)],[e5(:,1) e5(:,3)],[e9(:,1) e9(:,3)],[e6(:,1) e6(:,3)],[c1(1,1) c1(1,3)],[c2(1,1) c2(1,3)],[c6(1,1) c6(1,3)],[c5(1,1) c5(1,3)]);
f2(:,1:3) = [lattice(:,5) c1(1,2)*ones(length(lattice),1) lattice(:,6)];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
for i=1:Neta3
    f3((i-1)*Neta2+1:i*Neta2,1:3) = [e2(:,1)+e6(i,1)-c2(1,1) e2(:,2) zfunc(R,phi(i),0)*ones(length(e2),1)];
end

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% deltaz = (c6(1,3)-c2(1,3))/(Neta3-1);
% for i=1:Neta3
%     f3((i-1)*Neta2+1:i*Neta2,1:3) = [c2(1,1)*ones(length(ys),1) ys (c2(1,3)+deltaz*(i-1))*ones(length(ys),1)];
% end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e3(:,1) e3(:,3)],[e8(:,1) e8(:,3)],[e11(:,1) e11(:,3)],[e7(:,1) e7(:,3)],[c4(1,1) c4(1,3)],[c3(1,1) c3(1,3)],[c7(1,1) c7(1,3)],[c8(1,1) c8(1,3)]);
f4(:,1:3) = [lattice(:,5) c3(1,2)*ones(length(lattice),1) lattice(:,6)];

% ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
% deltaz = (c5(1,3)-c1(1,3))/(Neta3-1);
% for i=1:Neta3
%     f5((i-1)*Neta2+1:i*Neta2,1:3) = [c1(1,1)*ones(length(ys),1) ys (c1(1,3)+deltaz*(i-1))*ones(length(ys),1)];
% end

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
for i=1:Neta3
    f5((i-1)*Neta2+1:i*Neta2,1:3) = [e4(:,1)+e5(i,1)-c1(1,1) e4(:,2) zfunc(R,phi(i),0)*ones(length(e4),1)];
end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(logfullfile,Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e9,e12,e11,e10,c5,c6,c7,c8);
f6(:,1:3) = [lattice(:,5) lattice(:,6) c5(1,3)*ones(length(lattice),1)];

fig0 = figure();
plot3(c1(1,1),c1(1,2),c1(1,3),'kx')
hold on
plot3(c2(1,1),c2(1,2),c2(1,3),'kx')
hold on
plot3(c3(1,1),c3(1,2),c3(1,3),'kx')
hold on
plot3(c4(1,1),c4(1,2),c4(1,3),'kx')
hold on
plot3(c5(1,1),c5(1,2),c5(1,3),'kx')
hold on
plot3(c6(1,1),c6(1,2),c6(1,3),'kx')
hold on
plot3(c7(1,1),c7(1,2),c7(1,3),'kx')
hold on
plot3(c8(1,1),c8(1,2),c8(1,3),'kx')
hold on
plot3(e1(:,1),e1(:,2),e1(:,3),'r*')
hold on
plot3(e2(:,1),e2(:,2),e2(:,3),'r*')
hold on
plot3(e3(:,1),e3(:,2),e3(:,3),'r*')
hold on
plot3(e4(:,1),e4(:,2),e4(:,3),'r*')
hold on
plot3(e5(:,1),e5(:,2),e5(:,3),'r*')
hold on
plot3(e6(:,1),e6(:,2),e6(:,3),'r*')
hold on
plot3(e7(:,1),e7(:,2),e7(:,3),'r*')
hold on
plot3(e8(:,1),e8(:,2),e8(:,3),'r*')
hold on
plot3(e9(:,1),e9(:,2),e9(:,3),'r*')
hold on
plot3(e10(:,1),e10(:,2),e10(:,3),'r*')
hold on
plot3(e11(:,1),e11(:,2),e11(:,3),'r*')
hold on
plot3(e12(:,1),e12(:,2),e12(:,3),'r*')
hold on
plot3(f1(:,1),f1(:,2),f1(:,3),'b.')
hold on
plot3(f2(:,1),f2(:,2),f2(:,3),'b.')
hold on
plot3(f3(:,1),f3(:,2),f3(:,3),'b.')
hold on
plot3(f4(:,1),f4(:,2),f4(:,3),'b.')
hold on
plot3(f5(:,1),f5(:,2),f5(:,3),'b.')
hold on
plot3(f6(:,1),f6(:,2),f6(:,3),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

geom = generatelattice3D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,ymin,ymax,Ny,zmin,zmax,Nz);

geom = transfiniteinterpolation3D(Neta1,Neta2,Neta3,geom(:,1:3),eta1min,eta1max,eta2min,eta2max,eta3min,eta3max,1,f1,f2,f3,f4,f5,f6,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,c1,c2,c3,c4,c5,c6,c7,c8);

fig1 = figure();
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

fig2 = figure();
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
axis equal

fig3 = figure();
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(fig3,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

Neta1 = 30;
eta1min = 0;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

Neta2 = 25;
eta2min = 0;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

Neta3 = 35;
eta3min = 0;
deltaeta3 = 1;
eta3max = eta3min + (Neta3-1)*deltaeta3;

Nx = Neta1;
xmin = eta1min;
xmax = eta1max;

Ny = Neta2;
ymin = eta2min;
ymax = eta2max;

Nz = Neta3;
zmin = eta3min;
zmax = eta3max;

c1 = zeros(1,3);
c2 = zeros(1,3);
c3 = zeros(1,3);
c4 = zeros(1,3);
c5 = zeros(1,3);
c6 = zeros(1,3);
c7 = zeros(1,3);
c8 = zeros(1,3);

e1  = zeros(Neta1,3);
e2  = zeros(Neta2,3);
e3  = zeros(Neta1,3);
e4  = zeros(Neta2,3);
e5  = zeros(Neta3,3);
e6  = zeros(Neta3,3);
e7  = zeros(Neta3,3);
e8  = zeros(Neta3,3);
e9  = zeros(Neta1,3);
e10 = zeros(Neta2,3);
e11 = zeros(Neta1,3);
e12 = zeros(Neta2,3);

f1 = zeros(Neta1*Neta2,3);
f2 = zeros(Neta1*Neta3,3);
f3 = zeros(Neta2*Neta3,3);
f4 = zeros(Neta1*Neta3,3);
f5 = zeros(Neta2*Neta3,3);
f6 = zeros(Neta1*Neta2,3);

R = 20;

xfunc = @(r,phi,theta) r.*sin(phi).*cos(theta);
yfunc = @(r,phi,theta) r.*sin(phi).*sin(theta);
zfunc = @(r,phi,theta) r.*cos(phi);

c1(1,:) = [xfunc(R,0.75*pi,5*pi/4) yfunc(R,0.75*pi,5*pi/4) zfunc(R,0.75*pi,5*pi/4)];
c2(1,:) = [xfunc(R,0.75*pi,7*pi/4) yfunc(R,0.75*pi,7*pi/4) zfunc(R,0.75*pi,7*pi/4)];
c3(1,:) = [xfunc(R,0.75*pi,pi/4)   yfunc(R,0.75*pi,pi/4)   zfunc(R,0.75*pi,pi/4)];
c4(1,:) = [xfunc(R,0.75*pi,3*pi/4) yfunc(R,0.75*pi,3*pi/4) zfunc(R,0.75*pi,3*pi/4)];
c5(1,:) = [xfunc(R,0.25*pi,5*pi/4) yfunc(R,0.25*pi,5*pi/4) zfunc(R,0.25*pi,5*pi/4)];
c6(1,:) = [xfunc(R,0.25*pi,7*pi/4) yfunc(R,0.25*pi,7*pi/4) zfunc(R,0.25*pi,7*pi/4)];
c7(1,:) = [xfunc(R,0.25*pi,pi/4)   yfunc(R,0.25*pi,pi/4)   zfunc(R,0.25*pi,pi/4)];
c8(1,:) = [xfunc(R,0.25*pi,3*pi/4) yfunc(R,0.25*pi,3*pi/4) zfunc(R,0.25*pi,3*pi/4)];

rad = 0.5*abs(max(e6(:,1))-max(e5(:,1)));
x0 = 0.25*(c1(1,1)+c2(1,1)+c5(1,1)+c6(1,1));
z0 = 0.25*(c1(1,2)+c2(1,2)+c5(1,2)+c6(1,2));
theta = (5*pi/4:0.5*pi/(Neta1-1):7*pi/4)';
e1(:,1:3) = [x0+rad.*cos(theta) c1(1,2)*ones(Neta1,1) z0+rad*sin(theta)];

% xs = (c1(1,1):(c2(1,1)-c1(1,1))/(Neta1-1):c2(1,1))';
% e1(:,1:3) = [xs c1(1,2)*ones(length(xs),1) c1(1,3)*ones(length(xs),1)];

phi = 0.75*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e2(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% e2(:,1:3) = [c2(1,1)*ones(length(ys),1) ys c2(1,3)*ones(length(ys),1)];

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e3(:,1:3) = [xs c4(1,2)*ones(length(xs),1) c4(1,3)*ones(length(xs),1)];

phi = 0.75*pi;
theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e4(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
% e4(:,1:3) = [c1(1,1)*ones(length(ys),1) ys c1(1,3)*ones(length(ys),1)];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = pi;
e5(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c1(1,1) c1(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c1(1,3):(c5(1,3)-c1(1,3))/(Neta3-1):c5(1,3))';
% e5(:,1:3) = [c1(1,1)*ones(length(zs),1) c1(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e6(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c2(1,1) c2(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c2(1,3):(c6(1,3)-c2(1,3))/(Neta3-1):c6(1,3))';
% e6(:,1:3) = [c2(1,1)*ones(length(zs),1) c2(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = 0;
e7(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c3(1,1) c3(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c3(1,3):(c7(1,3)-c3(1,3))/(Neta3-1):c7(1,3))';
% e7(:,1:3) = [c3(1,1)*ones(length(zs),1) c3(1,2)*ones(length(zs),1) zs];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
theta = pi;
e8(:,1:3) = [xfunc(R,phi,theta)-R*sin(phi(1))*cos(theta)+c4(1,1) c4(1,2)*ones(length(phi),1) zfunc(R,phi,theta)];

% zs = (c4(1,3):(c8(1,3)-c4(1,3))/(Neta3-1):c8(1,3))';
% e8(:,1:3) = [c4(1,1)*ones(length(zs),1) c4(1,2)*ones(length(zs),1) zs];

xs = (c5(1,1):(c6(1,1)-c5(1,1))/(Neta1-1):c6(1,1))';
e9(:,1:3) = [xs c5(1,2)*ones(length(xs),1) c5(1,3)*ones(length(xs),1)];

phi = 0.25*pi;
theta = (-pi/4:0.5*pi/(Neta2-1):pi/4)';
e10(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c6(1,2):(c7(1,2)-c6(1,2))/(Neta2-1):c7(1,2))';
% e10(:,1:3) = [c6(1,1)*ones(length(ys),1) ys c6(1,3)*ones(length(ys),1)];

xs = (c8(1,1):(c7(1,1)-c8(1,1))/(Neta1-1):c7(1,1))';
e11(:,1:3) = [xs c8(1,2)*ones(length(xs),1) c8(1,3)*ones(length(xs),1)];

phi = 0.25*pi;
theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e12(:,1:3) = [xfunc(R,phi,theta) yfunc(R,phi,theta) zfunc(R,phi,theta).*ones(length(theta),1)];

% ys = (c5(1,2):(c8(1,2)-c5(1,2))/(Neta2-1):c8(1,2))';
% e12(:,1:3) = [c5(1,1)*ones(length(ys),1) ys c5(1,3)*ones(length(ys),1)];

p1 = 0.999;
p2 = 0.001;


lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e4,e3,e2,c1,c2,c3,c4);
f1(:,1:3) = [lattice(:,5) lattice(:,6) c1(1,3)*ones(length(lattice),1)];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e1(:,1) e1(:,3)],[e5(:,1) e5(:,3)],[e9(:,1) e9(:,3)],[e6(:,1) e6(:,3)],[c1(1,1) c1(1,3)],[c2(1,1) c2(1,3)],[c6(1,1) c6(1,3)],[c5(1,1) c5(1,3)]);
f2(:,1:3) = [lattice(:,5) c1(1,2)*ones(length(lattice),1) lattice(:,6)];

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
for i=1:Neta3
    f3((i-1)*Neta2+1:i*Neta2,1:3) = [e2(:,1)+e6(i,1)-c2(1,1) e2(:,2) zfunc(R,phi(i),0)*ones(length(e2),1)];
end

% ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
% deltaz = (c6(1,3)-c2(1,3))/(Neta3-1);
% for i=1:Neta3
%     f3((i-1)*Neta2+1:i*Neta2,1:3) = [c2(1,1)*ones(length(ys),1) ys (c2(1,3)+deltaz*(i-1))*ones(length(ys),1)];
% end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,zmin,zmax,Nz);
lattice = transfiniteinterpolation2D(Neta1*Neta3,lattice(:,1:2),eta1min,eta1max,eta3min,eta3max,Neta1,Neta3,1,[e3(:,1) e3(:,3)],[e8(:,1) e8(:,3)],[e11(:,1) e11(:,3)],[e7(:,1) e7(:,3)],[c4(1,1) c4(1,3)],[c3(1,1) c3(1,3)],[c7(1,1) c7(1,3)],[c8(1,1) c8(1,3)]);
f4(:,1:3) = [lattice(:,5) c3(1,2)*ones(length(lattice),1) lattice(:,6)];

% ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
% deltaz = (c5(1,3)-c1(1,3))/(Neta3-1);
% for i=1:Neta3
%     f5((i-1)*Neta2+1:i*Neta2,1:3) = [c1(1,1)*ones(length(ys),1) ys (c1(1,3)+deltaz*(i-1))*ones(length(ys),1)];
% end

phi = (0.75*pi:-0.5*pi/(Neta3-1):0.25*pi)';
for i=1:Neta3
    f5((i-1)*Neta2+1:i*Neta2,1:3) = [e4(:,1)+e5(i,1)-c1(1,1) e4(:,2) zfunc(R,phi(i),0)*ones(length(e4),1)];
end

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
lattice = transfiniteinterpolation2D(Neta1*Neta2,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e9,e12,e11,e10,c5,c6,c7,c8);
f6(:,1:3) = [lattice(:,5) lattice(:,6) c5(1,3)*ones(length(lattice),1)];

fig0 = figure();
plot3(c1(1,1),c1(1,2),c1(1,3),'kx')
hold on
plot3(c2(1,1),c2(1,2),c2(1,3),'kx')
hold on
plot3(c3(1,1),c3(1,2),c3(1,3),'kx')
hold on
plot3(c4(1,1),c4(1,2),c4(1,3),'kx')
hold on
plot3(c5(1,1),c5(1,2),c5(1,3),'kx')
hold on
plot3(c6(1,1),c6(1,2),c6(1,3),'kx')
hold on
plot3(c7(1,1),c7(1,2),c7(1,3),'kx')
hold on
plot3(c8(1,1),c8(1,2),c8(1,3),'kx')
hold on
plot3(e1(:,1),e1(:,2),e1(:,3),'r*')
hold on
plot3(e2(:,1),e2(:,2),e2(:,3),'r*')
hold on
plot3(e3(:,1),e3(:,2),e3(:,3),'r*')
hold on
plot3(e4(:,1),e4(:,2),e4(:,3),'r*')
hold on
plot3(e5(:,1),e5(:,2),e5(:,3),'r*')
hold on
plot3(e6(:,1),e6(:,2),e6(:,3),'r*')
hold on
plot3(e7(:,1),e7(:,2),e7(:,3),'r*')
hold on
plot3(e8(:,1),e8(:,2),e8(:,3),'r*')
hold on
plot3(e9(:,1),e9(:,2),e9(:,3),'r*')
hold on
plot3(e10(:,1),e10(:,2),e10(:,3),'r*')
hold on
plot3(e11(:,1),e11(:,2),e11(:,3),'r*')
hold on
plot3(e12(:,1),e12(:,2),e12(:,3),'r*')
hold on
% plot3(f1(:,1),f1(:,2),f1(:,3),'b.')
% hold on
% plot3(f2(:,1),f2(:,2),f2(:,3),'b.')
% hold on
% plot3(f3(:,1),f3(:,2),f3(:,3),'b.')
% hold on
% plot3(f4(:,1),f4(:,2),f4(:,3),'b.')
% hold on
% plot3(f5(:,1),f5(:,2),f5(:,3),'b.')
% hold on
% plot3(f6(:,1),f6(:,2),f6(:,3),'b.')
% hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

geom = generatelattice3D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,ymin,ymax,Ny,zmin,zmax,Nz);

geom = transfiniteinterpolation3D(Neta1,Neta2,Neta3,geom(:,1:3),eta1min,eta1max,eta2min,eta2max,eta3min,eta3max,1,f1,f2,f3,f4,f5,f6,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,c1,c2,c3,c4,c5,c6,c7,c8);

fig1 = figure();
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Sphere')
axis equal

fig2 = figure();
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
axis equal

fig3 = figure();
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(fig3,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(geom(:,1),geom(:,2),geom(:,3),'b.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(geom(:,7),geom(:,8),geom(:,9),'b.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
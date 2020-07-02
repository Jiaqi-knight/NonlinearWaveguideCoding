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

x1 = 7;
y1 = 0;
R1 = 12;

x2 = -9;
y2 = 0;
R2 = 17;

R3 = 3;

d = sqrt((x1-x2)^2+(y1-y2)^2);

delta = 0.5*d+0.5*(R2^2-R1^2)/d;
xa = x2 + delta;
xb = xa;

ya = y2 + sqrt(R2^2-delta^2);
yb = y2 - sqrt(R2^2-delta^2);

yc = ya;
err = 1;
it = 0;
tol = 10^-14;
itmax = 100;
while err>=tol && it<itmax
    f = sqrt((R1+R3)^2-yc^2) + sqrt((R2+R3)^2-yc^2) -d;
    df = -yc*(1/sqrt((R1+R3)^2-yc^2) + 1/sqrt((R2+R3)^2-yc^2));
    yc = yc - f/df;
    it = it + 1;
    err = abs(sqrt((R1+R3)^2-yc^2) + sqrt((R2+R3)^2-yc^2) -d);
end

yc = real(yc);
yd = -yc;

xc = x2 + sqrt((R2+R3)^2-yc^2);
xd = xc;

beta = asin(yc/(R1+R3));
gamma = asin(yc/(R2+R3));

xq = x1 - R1*cos(beta);
yq = y1 + R1*sin(beta);

xp = x2 + R2*cos(gamma);
yp = y2 + R2*sin(gamma);

base1 = sqrt(R2^2-ya^2);
base2 = sqrt(R1^2-ya^2);

phi = atan(ya/base1);
psi = atan(ya/base2);

eta1min = 0;
Neta1_1 = 20;
Neta1_2 = 5;
Neta1_3 = 5;
Neta1 = Neta1_1 + Neta1_2 + Neta1_3;
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

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);

c1(1,1) = x2;
c1(1,2) = y2-R2;

c2(1,1) = x1;
c2(1,2) = y1-R1;

c3(1,1) = x1;
c3(1,2) = y1+R1;

c4(1,1) = x2;
c4(1,2) = y2+R2;

theta = (-0.5*pi:pi/(Neta2-1):0.5*pi)';
e4(:,1) = x1 + R1*cos(theta);
e4(:,2) = y1 + R1*sin(theta);

theta = (1.5*pi:-pi/(Neta2-1):0.5*pi)';
e2(:,1) = x2 + R2*cos(theta);
e2(:,2) = y2 + R2*sin(theta);

gammacomp = 0.5*pi - gamma;
betacomp = 0.5*pi - beta;

theta1 = (0.5*pi:(gamma-0.5*pi)/(Neta1_1-1):gamma)';
theta2 = (pi+gammacomp:((1.5*pi+betacomp)-(pi+gammacomp))/(Neta1_2-1):1.5*pi+betacomp)';
theta3 = (pi-beta:-(0.5*pi-beta)/(Neta1_3-1):0.5*pi)';
e3(:,1) = [x2+R2.*cos(theta1); xc+R3.*cos(theta2); x1+R1.*cos(theta3)];
e3(:,2) = [y2+R2.*sin(theta1); yc+R3.*sin(theta2); y1+R1.*sin(theta3)];

e1(:,1) = e3(:,1);
e1(:,2) = -e3(:,2);

f0 = figure();
plot(x1,y1,'g*','LineWidth',2)
hold on
plot(x2,y2,'g*','LineWidth',2)
hold on
plot([x1;xc],[y1;yc],'b--','LineWidth',2)
hold on
plot([x2;xc],[y2;yc],'b--','LineWidth',2)
hold on
plot([x1;xq],[y1;yq],'b--','LineWidth',2)
hold on
plot([x2;xp],[y2;yp],'b--','LineWidth',2)
hold on
plot([x2;x1],[y2;y1],'b--','LineWidth',2)
hold on
plot([xa;xa],[yc+5;yd-5],'k--','LineWidth',2)
hold on
plot([x2-R2-5;x1+R1+5],[y2;y1],'k--','LineWidth',2)
hold on
plot(xa,ya,'g*','LineWidth',2)
hold on
plot(xb,yb,'g*','LineWidth',2)
hold on
plot(xc,yc,'g*','LineWidth',2)
hold on
plot(xd,yd,'g*','LineWidth',2)
hold on
plot(xp,yp,'g*','LineWidth',2)
hold on
plot(xq,yq,'g*','LineWidth',2)
hold on
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

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

lattice = transfiniteinterpolation2D(logfullfile,N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,Nx,Ny);

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

figure()
plot(lattice(:,3),lattice(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

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
%plot(lattice(indicesE2,3),lattice(indicesE2,4),'.')
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

figure();
plot(lattice(indicesC1,1),lattice(indicesC1,2),'.')
hold on
plot(lattice(indicesC2,1),lattice(indicesC2,2),'.')
hold on
plot(lattice(indicesC3,1),lattice(indicesC3,2),'.')
hold on
plot(lattice(indicesC4,1),lattice(indicesC4,2),'.')
hold on
plot(lattice(indicesE1,1),lattice(indicesE1,2),'.')
hold on
%plot(lattice(indicesE2,1),lattice(indicesE2,2),'.')
hold on
plot(lattice(indicesE3,1),lattice(indicesE3,2),'.')
hold on
plot(lattice(indicesE4,1),lattice(indicesE4,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - transfinite interpolation')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periodicity = 0;
deltaq = [deltaeta1 deltaeta2];
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(logfullfile,N,Nx,0,0,0,0,0,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
                                                                                                  
itmax = 10;
tol = 10^-4;

lattice10 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice50 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice100 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice1000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice5000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice10000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice50000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

lattice100000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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

itmax = 500000;
tol = 10^-4;

lattice500000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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
plot(lattice500000(:,1),lattice500000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice500000(:,5),lattice500000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 500000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 1000000;
tol = 10^-4;

lattice1000000 = ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol); 

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
plot(lattice1000000(:,1),lattice1000000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice1000000(:,5),lattice1000000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1000000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periodicity = 0;
deltaq = [deltaeta1 deltaeta2];
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(logfullfile,N,Nx,0,0,0,0,0,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

itmax = 10;
tol = 10^-4;

lattice10 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice50 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice100 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 
             %sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,flagperiodicity,periodicity,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,spyflag)
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

lattice1000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice5000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice10000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice50000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

lattice100000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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

itmax = 500000;
tol = 10^-4;

lattice500000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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
plot(lattice500000(:,1),lattice500000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice500000(:,5),lattice500000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 500000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])

itmax = 1000000;
tol = 10^-4;

lattice1000000 = sparseellipticgridgen2D(logfullfile,Nx,N,lattice,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

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
plot(lattice1000000(:,1),lattice1000000(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice1000000(:,5),lattice1000000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1000000 it')

axis([fcomp fphys],'equal')
axis(fcomp,[0 100 0 100])
%axis(fphys,[0 11 -5 11])
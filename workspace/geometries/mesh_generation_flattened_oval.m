%clear all
%close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = 7;
y1 = 0;
R1 = 17;

x2 = -9;
y2 = 0;
R2 = 17;

d = sqrt((x1-x2)^2+(y1-y2)^2);

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

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

c1(1,1:2) = [x2 y2-R2];
c2(1,1:2) = [x1 y1-R1];
c3(1,1:2) = [x1 y1+R1];
c4(1,1:2) = [x2 y2+R2];

xs = (c1(1):(c2(1)-c1(1))/(Neta1-1):c2(1))';
me1 = (c2(2)-c1(2))/(c2(1)-c1(1));
e1(:,1) = xs;
e1(:,2) = me1.*(xs-c1(1))+c1(2);

theta = (-0.5*pi:pi/(Neta2-1):0.5*pi)';
e4(:,1) = x1+R1*cos(theta);
e4(:,2) = y1+R1*sin(theta);

xs = (c4(1):(c3(1)-c4(1))/(Neta1-1):c3(1))';
me3 = (c3(2)-c4(2))/(c3(1)-c4(1));
e3(:,1) = xs;
e3(:,2) = me3.*(xs-c4(1))+c4(2);

theta = (1.5*pi:-pi/(Neta2-1):0.5*pi)';
e2(:,1) = x2+R2*cos(theta);
e2(:,2) = y2+R2*sin(theta);


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

lattice = transfiniteinterpolation2D(N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

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


figure();
plot(lattice(:,3),lattice(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50;

x1 = 1;
x2 = 5;
x3_4 = 12;

c = 10;

[airfoil,yc] = NACAxxxx(N,x1,x2,x3_4,c);

figure();
plot(airfoil(:,1),airfoil(:,2),'*')
hold on
plot(airfoil(1:length(yc),1),yc,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))

figure();
plot(airfoil(1:length(yc),1),airfoil(1:length(yc),2),'*')
hold on
plot(airfoil(1:length(yc),1),yc,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))

figure();
plot(airfoil(length(yc)+1:end,1),airfoil(length(yc)+1:end,2),'*')
hold on
plot(airfoil(1:length(yc),1),yc,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))

Nleading = 50;
Ntrailing = 125;

xmin = min(airfoil(:,1))-3*c;
xmax = max(airfoil(:,1))+6*c;

xleadinglower = (xmin:(min(airfoil(:,1))-xmin)/Nleading:min(airfoil(:,1)))';
xleadingupper = xleadinglower(1:end-1,1);

xtrailinglower = (max(airfoil(:,1)):(xmax-max(airfoil(:,1)))/Ntrailing:xmax)';
xtrailingupper = xtrailinglower(2:end,1);

upper = [xleadingupper yc(1)*ones(length(xleadingupper),1);airfoil(1:length(yc),1) airfoil(1:length(yc),2);xtrailingupper yc(end)*ones(length(xtrailingupper),1)];

lower = [xleadinglower airfoil(length(yc)+2,2)*ones(length(xleadinglower),1);airfoil(length(yc)+2:end-1,1) airfoil(length(yc)+2:end-1,2);xtrailinglower airfoil(end-1,2)*ones(length(xtrailinglower),1)];

figure();
plot(upper(:,1),upper(:,2),'*')
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))

figure();
plot(lower(:,1),lower(:,2),'*')
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))

h = 10;

eta1min = 0;
Neta1 = Nleading + N + Ntrailing;
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

c1(1,1) = upper(1,1);
c1(1,2) = upper(1,2);

c2(1,1) = upper(end,1);
c2(1,2) = upper(end,2);

c3(1,1) = upper(end,1);
c3(1,2) = max(upper(1,2),upper(end,2))+h;

c4(1,1) = upper(1,1);
c4(1,2) = max(upper(1,2),upper(end,2))+h;

e1(:,1) = upper(:,1);
e1(:,2) = upper(:,2);

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
e2(:,1) = upper(1,1)*ones(length(ys),1);
e2(:,2) = ys;

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e3(:,1) = xs;
e3(:,2) = c4(1,2)*ones(length(xs),1);

ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
e4(:,1) = upper(end,1)*ones(length(ys),1);
e4(:,2) = ys;

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

latticeupper = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

latticeupper = transfiniteinterpolation2D(logfullfile,N,latticeupper(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulkupper,indicesinternalbulkupper,indicesE1upper,indicesE2upper,indicesE3upper,indicesE4upper,indicesexternalE1upper,indicesexternalE2upper,indicesexternalE3upper,indicesexternalE4upper,indicesinternalE1upper,indicesinternalE2upper,indicesinternalE3upper,indicesinternalE4upper,indicesC1upper,indicesC2upper,indicesC3upper,indicesC4upper,indicesinternalC1upper,indicesinternalC2upper,indicesinternalC3upper,indicesinternalC4upper] = getindices2D(logfullfile,Nx,Ny);

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticeupper(:,3),latticeupper(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

N = 50;

eta1min = 0;
Neta1 = Nleading + N + Ntrailing;
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

c1(1,1) = lower(1,1);
c1(1,2) = min(lower(1,2),lower(end,2))-h;

c2(1,1) = lower(end,1);
c2(1,2) = min(lower(1,2),lower(end,2))-h;

c3(1,1) = lower(end,1);
c3(1,2) = lower(end,2);

c4(1,1) = lower(1,1);
c4(1,2) = lower(1,2);

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e1(:,1) = xs;
e1(:,2) = c1(1,2)*ones(length(xs),1);

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
e2(:,1) = lower(1,1)*ones(length(ys),1);
e2(:,2) = ys;

e3(:,1) = lower(:,1);
e3(:,2) = lower(:,2);

ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
e4(:,1) = lower(end,1)*ones(length(ys),1);
e4(:,2) = ys;

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

latticelower = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

latticelower = transfiniteinterpolation2D(logfullfile,N,latticelower(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulklower,indicesinternalbulklower,indicesE1lower,indicesE2lower,indicesE3lower,indicesE4lower,indicesexternalE1lower,indicesexternalE2lower,indicesexternalE3lower,indicesexternalE4lower,indicesinternalE1lower,indicesinternalE2lower,indicesinternalE3lower,indicesinternalE4lower,indicesC1lower,indicesC2lower,indicesC3lower,indicesC4lower,indicesinternalC1lower,indicesinternalC2lower,indicesinternalC3lower,indicesinternalC4lower] = getindices2D(logfullfile,Nx,Ny);

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticelower(:,3),latticelower(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticeupper(:,3),latticeupper(:,4),'.')
hold on
plot(latticelower(:,3),latticelower(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

lattice = [latticelower;latticeupper(:,1) deltaeta2*Neta2+latticeupper(:,2) latticeupper(:,3:end)];

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
plot(airfoil(:,1),airfoil(:,2),'*r')
hold on
plot(airfoil(1:length(yc),1),yc,'k--','LineWidth',2)
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

Neta2new = 2*Neta2;

lowerstartindex = length(xleadinglower)+(Neta2-1)*Neta1;
upperstartindex = length(xleadingupper)+Neta2*Neta1;

lowernum = length(yc)-2;
uppernum = length(yc);

boundary = zeros(2*length(yc)-2,1);
for i=1:lowernum
    boundary(i,1) = lowerstartindex + i;
end
for i=1:uppernum
    boundary(lowernum+i,1) = upperstartindex + i;
end

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice(boundary,1),lattice(boundary,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice(boundary,3),lattice(boundary,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(lattice(boundary,1),lattice(boundary,2),'.r')
hold on
for i=1:size(lattice,1)
   if ~any(boundary==i)
       plot(lattice(i,1),lattice(i,2),'.b')
       hold on
   end
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(lattice(boundary,3),lattice(boundary,4),'.r')
hold on
for i=1:size(lattice,1)
   if ~any(boundary==i)
       plot(lattice(i,3),lattice(i,4),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% HORIZONTAL TAIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nht = 25;

x1ht = 0;
x2ht = 0;
x3_4ht = 08;

cht = 4;

hht = 10;

[horztail,ycht] = NACAxxxx(Nht,x1ht,x2ht,x3_4ht,cht);

figure();
plot(horztail(:,1),horztail(:,2),'*')
hold on
%plot(horztail(1:length(ycht),1),ycht,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1ht),num2str(x2ht),num2str(x3_4ht),' profile'))
axis equal

figure();
plot(horztail(1:length(ycht),1),horztail(1:length(ycht),2),'*')
hold on
plot(horztail(1:length(ycht),1),ycht,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1ht),num2str(x2ht),num2str(x3_4ht),' profile'))
axis equal

figure();
plot(horztail(length(ycht)+1:end,1),horztail(length(ycht)+1:end,2),'*')
hold on
plot(horztail(1:length(ycht),1),ycht,'k--','LineWidth',2)
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1ht),num2str(x2ht),num2str(x3_4ht),' profile'))
axis equal

x0 = 20;
y0 = 15;

Nleadinght = 150;
Ntrailinght = 50;

xminht = min(airfoil(:,1))-3*c;
xmaxht = max(airfoil(:,1))+6*c;

xleadinglowerht = (xminht:(x0+min(horztail(:,1))-xminht)/Nleadinght:x0+min(horztail(:,1)))';
xleadingupperht = xleadinglowerht(1:end-1,1);

xtrailinglowerht = (x0+max(horztail(:,1)):(xmaxht-(x0+max(horztail(:,1))))/Ntrailinght:xmaxht)';
xtrailingupperht = xtrailinglowerht(2:end,1);

upperht = [xleadingupperht y0+ycht(1)*ones(length(xleadingupperht),1);x0+horztail(1:length(ycht),1) y0+horztail(1:length(ycht),2);xtrailingupperht y0+ycht(end)*ones(length(xtrailingupperht),1)];

lowerht = [xleadinglowerht y0+horztail(length(ycht)+2,2)*ones(length(xleadinglowerht),1);x0+horztail(length(ycht)+2:end-1,1) y0+horztail(length(ycht)+2:end-1,2);xtrailinglowerht y0+horztail(end-1,2)*ones(length(xtrailinglowerht),1)];

figure();
plot(upperht(:,1),upperht(:,2),'*')
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1ht),num2str(x2ht),num2str(x3_4ht),' profile'))

figure();
plot(lowerht(:,1),lowerht(:,2),'*')
hold on
grid on
axis equal
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1ht),num2str(x2ht),num2str(x3_4ht),' profile'))

eta1min = 0;
Neta1ht = Nleadinght + Nht + Ntrailinght;
deltaeta1 = 1;
eta1max = eta1min + Neta1ht*deltaeta1;
eta2min = 0;
Neta2ht = 60;
deltaeta2 = 1;
eta2max = eta2min + Neta2ht*deltaeta2;

N = Neta1ht*Neta2ht;

xmin = 0;
xmax = 5;
Nx = Neta1ht;
ymin = -2.5;
ymax = 3;
Ny = Neta2ht;

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1ht,6);
e2 = zeros(Neta2ht,6);
e3 = zeros(Neta1ht,6);
e4 = zeros(Neta2ht,6);

c1(1,1) = upperht(1,1);
c1(1,2) = upperht(1,2);

c2(1,1) = upperht(end,1);
c2(1,2) = upperht(end,2);

c3(1,1) = upperht(end,1);
c3(1,2) = max(upperht(1,2),upperht(end,2))+hht;

c4(1,1) = upperht(1,1);
c4(1,2) = max(upperht(1,2),upperht(end,2))+hht;

e1(:,1) = upperht(:,1);
e1(:,2) = upperht(:,2);

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2ht-1):c4(1,2))';
e2(:,1) = upperht(1,1)*ones(length(ys),1);
e2(:,2) = ys;

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1ht-1):c3(1,1))';
e3(:,1) = xs;
e3(:,2) = c4(1,2)*ones(length(xs),1);

ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2ht-1):c3(1,2))';
e4(:,1) = upperht(end,1)*ones(length(ys),1);
e4(:,2) = ys;

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

latticeupperht = generatelattice2D(eta1min,Neta1ht,deltaeta1,eta2min,Neta2ht,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

latticeupperht = transfiniteinterpolation2D(N,latticeupperht(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1ht,Neta2ht,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulkupper,indicesinternalbulkupper,indicesE1upper,indicesE2upper,indicesE3upper,indicesE4upper,indicesexternalE1upper,indicesexternalE2upper,indicesexternalE3upper,indicesexternalE4upper,indicesinternalE1upper,indicesinternalE2upper,indicesinternalE3upper,indicesinternalE4upper,indicesC1upper,indicesC2upper,indicesC3upper,indicesC4upper,indicesinternalC1upper,indicesinternalC2upper,indicesinternalC3upper,indicesinternalC4upper] = getindices2D(Nx,Ny);

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticeupperht(:,3),latticeupperht(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal


eta1min = 0;
Neta1ht = Nleadinght + Nht + Ntrailinght;
deltaeta1 = 1;
eta1max = eta1min + Neta1ht*deltaeta1;
eta2min = 0;
Neta2ht = 60;
deltaeta2 = 1;
eta2max = eta2min + Neta2ht*deltaeta2;

N = Neta1ht*Neta2ht;

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

c1(1,1) = lowerht(1,1);
c1(1,2) = max(lattice(:,4))+0.5*(lattice(end,4)-lattice(end-Neta1,4));

c2(1,1) = lowerht(end,1);
c2(1,2) = max(lattice(:,4))+0.5*(lattice(end,4)-lattice(end-Neta1,4));

c3(1,1) = lowerht(end,1);
c3(1,2) = lowerht(end,2);

c4(1,1) = lowerht(1,1);
c4(1,2) = lowerht(1,2);

xs = (c4(1,1):(c3(1,1)-c4(1,1))/(Neta1-1):c3(1,1))';
e1(:,1) = xs;
e1(:,2) = c1(1,2)*ones(length(xs),1);

ys = (c1(1,2):(c4(1,2)-c1(1,2))/(Neta2-1):c4(1,2))';
e2(:,1) = lowerht(1,1)*ones(length(ys),1);
e2(:,2) = ys;

e3(:,1) = lowerht(:,1);
e3(:,2) = lowerht(:,2);

ys = (c2(1,2):(c3(1,2)-c2(1,2))/(Neta2-1):c3(1,2))';
e4(:,1) = lowerht(end,1)*ones(length(ys),1);
e4(:,2) = ys;

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

latticelowerht = generatelattice2D(eta1min,Neta1ht,deltaeta1,eta2min,Neta2ht,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

latticelowerht = transfiniteinterpolation2D(N,latticelowerht(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1ht,Neta2ht,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulklower,indicesinternalbulklower,indicesE1lower,indicesE2lower,indicesE3lower,indicesE4lower,indicesexternalE1lower,indicesexternalE2lower,indicesexternalE3lower,indicesexternalE4lower,indicesinternalE1lower,indicesinternalE2lower,indicesinternalE3lower,indicesinternalE4lower,indicesC1lower,indicesC2lower,indicesC3lower,indicesC4lower,indicesinternalC1lower,indicesinternalC2lower,indicesinternalC3lower,indicesinternalC4lower] = getindices2D(Nx,Ny);

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticelowerht(:,3),latticelowerht(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(latticeupperht(:,3),latticeupperht(:,4),'.')
hold on
plot(latticelowerht(:,3),latticelowerht(:,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

latticeht = [latticelowerht;latticeupperht(:,1) deltaeta2*Neta2+latticeupperht(:,2) latticeupperht(:,3:end)];

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
plot(latticeht(:,1),latticeht(:,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(latticeht(:,3),latticeht(:,4),'.')
hold on
plot(horztail(:,1),horztail(:,2),'*r')
hold on
plot(horztail(1:length(ycht),1),ycht,'k--','LineWidth',2)
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

Neta2newht = 2*Neta2ht;

lowerstartindexht = length(xleadinglowerht)+(Neta2ht-1)*Neta1ht;
upperstartindexht = length(xleadingupperht)+Neta2ht*Neta1ht;

lowernumht = length(ycht)-2;
uppernumht = length(ycht);

boundaryht = zeros(2*length(ycht)-2,1);
for i=1:lowernumht
    boundaryht(i,1) = lowerstartindexht + i;
end
for i=1:uppernumht
    boundaryht(lowernumht+i,1) = upperstartindexht + i;
end

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(latticeht(boundaryht,1),latticeht(boundaryht,2),'.')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(latticeht(boundaryht,3),latticeht(boundaryht,4),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(latticeht(boundaryht,1),latticeht(boundaryht,2),'.r')
hold on
for i=1:size(latticeht,1)
   if ~any(boundaryht==i)
       plot(latticeht(i,1),latticeht(i,2),'.b')
       hold on
   end
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(latticeht(boundaryht,3),latticeht(boundaryht,4),'.r')
hold on
for i=1:size(latticeht,1)
   if ~any(boundaryht==i)
       plot(latticeht(i,3),latticeht(i,4),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

latticeglobal = [lattice;latticeht(:,1) deltaeta2*Neta2new+latticeht(:,2) latticeht(:,3:end)];

boundaryhtglobal = size(lattice,1) + boundaryht;

boundaryglobal = [boundary;boundaryhtglobal];

f = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on

subplot(1,2,1);
subplot(1,2,2);

hsubfigs = get(f,'Children');

fcomp = hsubfigs(2);
fphys = hsubfigs(1);

subplot(fcomp);
plot(latticeglobal(boundaryglobal,1),latticeglobal(boundaryglobal,2),'.r')
hold on
for i=1:size(latticeglobal,1)
   if ~any(boundaryglobal==i)
       plot(latticeglobal(i,1),latticeglobal(i,2),'.b')
       hold on
   end
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')

subplot(fphys);
plot(latticeglobal(boundaryglobal,3),latticeglobal(boundaryglobal,4),'.r')
hold on
for i=1:size(latticeglobal,1)
   if ~any(boundaryglobal==i)
       plot(latticeglobal(i,3),latticeglobal(i,4),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis([fcomp fphys],'equal')

%close all

figure();
plot(latticeglobal(boundaryglobal,3),latticeglobal(boundaryglobal,4),'.r')
hold on
for i=1:size(latticeglobal,1)
   if ~any(boundaryglobal==i)
       plot(latticeglobal(i,3),latticeglobal(i,4),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

flagintbounds = 1;
indicesintbounds = boundaryglobal;
typeintbounds = [7*ones(lowernum,1);4;5*ones(uppernum-2,1);3;7*ones(lowernumht,1);4;5*ones(uppernumht-2,1);3];

clear latticeht latticelower latticelowerht latticeupper lattice latticeupperht upperht lower lowerht lowernum lowernumht lowerstartindex lowerstartindex lowerstartindexht
clear upper upperht uppernum uppernumht upperstartindex upperstartindex upperstartindexht xleadinglower xleadinglowerht xleadingupper xleadingupperht
clear xtrailinglower xtrailinglowerht xtrailingupper xtrailingupperht xs yc ycht ys airfoil horztail boundary boundaryht boundaryhtglobal e1 e2 e3 e4
clear c1 c2 c3 c4
clear indicesbulkupper indicesinternalbulkupper indicesE1upper indicesE2upper indicesE3upper indicesE4upper indicesexternalE1upper indicesexternalE2upper indicesexternalE3upper indicesexternalE4upper indicesinternalE1upper indicesinternalE2upper indicesinternalE3upper indicesinternalE4upper indicesC1upper indicesC2upper indicesC3upper indicesC4upper indicesinternalC1upper indicesinternalC2upper indicesinternalC3upper indicesinternalC4upper 
clear indicesbulklower indicesinternalbulklower indicesE1lower indicesE2lower indicesE3lower indicesE4lower indicesexternalE1lower indicesexternalE2lower indicesexternalE3lower indicesexternalE4lower indicesinternalE1lower indicesinternalE2lower indicesinternalE3lower indicesinternalE4lower indicesC1lower indicesC2lower indicesC3lower indicesC4lower indicesinternalC1lower indicesinternalC2lower indicesinternalC3lower indicesinternalC4lower 
clear eta2min eta2max f fcomp fphys h hht hsubfigs i N x0 x1 x1ht x2 x2ht x3_4 x3_4ht xmax xmaxht xmin xminht y0 ymax ymin scrsz c cht

flagperiodicity = 0;
periodicity = 0;
deltaq = [deltaeta1 deltaeta2];
[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Neta1,(Neta2new+2*Neta2));
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(Neta1*(Neta2new+2*Neta2),Neta1,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

clear indicesinternalbulk indicesexternalE1 indicesexternalE2 indicesexternalE3 indicesexternalE4 indicesinternalE1 indicesinternalE2 indicesinternalE3 indicesinternalE4 indicesinternalC1 indicesinternalC2 indicesinternalC3 indicesinternalC4
clear structuralneighbours shearneighbours bendneighbours

itmax = 1;
tol = 10^-4;

lattice1 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice1(:,5),lattice1(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1 it')
axis equal

figure();
plot(lattice1(boundaryglobal,5),lattice1(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice1,1)
   if ~any(boundaryglobal==i)
       plot(lattice1(i,5),lattice1(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice1

itmax = 2;
tol = 10^-4;

lattice2 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice2(:,5),lattice2(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 2 it')
axis equal

figure();
plot(lattice2(boundaryglobal,5),lattice2(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice2,1)
   if ~any(boundaryglobal==i)
       plot(lattice2(i,5),lattice2(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice2

itmax = 10;
tol = 10^-4;

lattice10 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice10(:,5),lattice10(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 10 it')
axis equal

figure();
plot(lattice10(boundaryglobal,5),lattice10(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice10,1)
   if ~any(boundaryglobal==i)
       plot(lattice10(i,5),lattice10(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice10

itmax = 100;
tol = 10^-4;

lattice100 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice100(:,5),lattice100(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 100 it')
axis equal

figure();
plot(lattice100(boundaryglobal,5),lattice100(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice100,1)
   if ~any(boundaryglobal==i)
       plot(lattice100(i,5),lattice100(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice100

itmax = 1000;
tol = 10^-4;

lattice1000 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice1000(:,5),lattice1000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1000 it')
axis equal

figure();
plot(lattice1000(boundaryglobal,5),lattice1000(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice1000,1)
   if ~any(boundaryglobal==i)
       plot(lattice1000(i,5),lattice1000(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice1000

itmax = 10000;
tol = 10^-4;

lattice10000 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice10000(:,5),lattice10000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 10000 it')
axis equal

figure();
plot(lattice10000(boundaryglobal,5),lattice10000(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice10000,1)
   if ~any(boundaryglobal==i)
       plot(lattice10000(i,5),lattice10000(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice10000

itmax = 100000;
tol = 10^-4;

lattice100000 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice100000(:,5),lattice100000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 100000 it')
axis equal

figure();
plot(lattice100000(boundaryglobal,5),lattice100000(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice100000,1)
   if ~any(boundaryglobal==i)
       plot(lattice100000(i,5),lattice100000(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal


clear lattice100000

itmax = 1000000;
tol = 10^-4;

lattice1000000 = sparseellipticgridgen2Dinternalboundaries(Neta1,Neta1*(Neta2new+2*Neta2),latticeglobal,deltaq,boundaryglobal,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 

figure();
plot(lattice1000000(:,5),lattice1000000(:,6),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space) - elliptic smoothing 1000000 it')
axis equal

figure();
plot(lattice1000000(boundaryglobal,5),lattice1000000(boundaryglobal,6),'.r')
hold on
for i=1:size(lattice1000000,1)
   if ~any(boundaryglobal==i)
       plot(lattice1000000(i,5),lattice1000000(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')

axis equal

clear lattice1000000 latticeglobal

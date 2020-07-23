clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DEFORMED MESHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta1min = 0;
Neta1 = 10;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

eta2min = 0;
Neta2 = 10;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

eta3min = 0;
Neta3 = 10;
deltaeta3 = 1;
eta3max = eta3min + (Neta3-1)*deltaeta3;

deltaq = [deltaeta1 deltaeta2];

Nx = Neta1;
xmin = 0;
xmax = 5;
deltax = (xmax-xmin)/(Nx-1);

Ny = Neta2;
ymin = -2.5;
ymax = 3;
deltay = (ymax-ymin)/(Ny-1);

N = Nx*Ny;

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

%point1 = Nx + floor(Nx/2) + randi(Nx*Ny-2*(Nx+floor(Nx/2)));
point1 = 65;
point2 = point1 + 1 + Nx;

scrsz = get(0,'ScreenSize');
f1 = figure();
for i=1:Nx*Ny
    if i~=point1 && i~=point2
        plot(lattice(i,5),lattice(i,6),'.b')
        hold on
    end
end
plot(lattice(point1,5),lattice(point1,6),'*r')
hold on
plot(lattice(point2,5),lattice(point2,6),'*g')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in physical space')
axis equal

scrsz = get(0,'ScreenSize');
f2 = figure();
for i=1:Nx*Ny
    if i~=point1 && i~=point2
        plot(lattice(i,1),lattice(i,2),'.b')
        hold on
    end
end
plot(lattice(point1,1),lattice(point1,2),'*r')
hold on
plot(lattice(point2,1),lattice(point2,2),'*g')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in computational space')
axis equal

temp = lattice(point2,5:6);
lattice(point2,5:6) = lattice(point1,5:6) + [rand()*0.5*deltax rand()*0.5*deltay];
lattice(point1,5:6) = temp + [rand()*0.5*deltax rand()*0.5*deltay];

scrsz = get(0,'ScreenSize');
f3 = figure();
for i=1:Nx*Ny
    if i~=point1 && i~=point2
        plot(lattice(i,5),lattice(i,6),'.b')
        hold on
    end
end
plot(lattice(point1,5),lattice(point1,6),'*r')
hold on
plot(lattice(point2,5),lattice(point2,6),'*g')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in physical space')
axis equal

scrsz = get(0,'ScreenSize');
f4 = figure();
for i=1:Nx*Ny
    if i~=point1 && i~=point2
        plot(lattice(i,1),lattice(i,2),'.b')
        hold on
    end
end
plot(lattice(point1,1),lattice(point1,2),'*r')
hold on
plot(lattice(point2,1),lattice(point2,2),'*g')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in computational space')
axis equal

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,Nx,Ny);

%% figure
figure();
for i=1:length(indicesbulk)
        plot(lattice(indicesbulk(i),1),lattice(indicesbulk(i),2),'.b')
        hold on
end
for i=1:length(indicesinternalbulk)
        plot(lattice(indicesinternalbulk(i),1),lattice(indicesinternalbulk(i),2),'.r')
        hold on
end
for i=1:length(indicesE1)
        plot(lattice(indicesE1(i),1),lattice(indicesE1(i),2),'*r')
        hold on
end
for i=1:length(indicesE2)
        plot(lattice(indicesE2(i),1),lattice(indicesE2(i),2),'*b')
        hold on
end
for i=1:length(indicesE3)
        plot(lattice(indicesE3(i),1),lattice(indicesE3(i),2),'*r')
        hold on
end
for i=1:length(indicesE4)
        plot(lattice(indicesE4(i),1),lattice(indicesE4(i),2),'*b')
        hold on
end
for i=1:length(indicesexternalE1)
        plot(lattice(indicesexternalE1(i),1),lattice(indicesexternalE1(i),2),'or')
        hold on
end
for i=1:length(indicesexternalE2)
        plot(lattice(indicesexternalE2(i),1),lattice(indicesexternalE2(i),2),'ob')
        hold on
end
for i=1:length(indicesexternalE3)
        plot(lattice(indicesexternalE3(i),1),lattice(indicesexternalE3(i),2),'or')
        hold on
end
for i=1:length(indicesexternalE4)
        plot(lattice(indicesexternalE4(i),1),lattice(indicesexternalE4(i),2),'ob')
        hold on
end
for i=1:length(indicesinternalE1)
        plot(lattice(indicesinternalE1(i),1),lattice(indicesinternalE1(i),2),'or')
        hold on
end
for i=1:length(indicesinternalE2)
        plot(lattice(indicesinternalE2(i),1),lattice(indicesinternalE2(i),2),'ob')
        hold on
end
for i=1:length(indicesinternalE3)
        plot(lattice(indicesinternalE3(i),1),lattice(indicesinternalE3(i),2),'or')
        hold on
end
for i=1:length(indicesinternalE4)
        plot(lattice(indicesinternalE4(i),1),lattice(indicesinternalE4(i),2),'ob')
        hold on
end
        plot(lattice(indicesC1,1),lattice(indicesC1,2),'^')
        plot(lattice(indicesC2,1),lattice(indicesC2,2),'^')
        plot(lattice(indicesC3,1),lattice(indicesC3,2),'^')
        plot(lattice(indicesC4,1),lattice(indicesC4,2),'^')
        
        plot(lattice(indicesinternalC1(1)-100,1),lattice(indicesinternalC1(1)-100,2),'Marker','pentagram')
        plot(lattice(indicesinternalC2(1)-100,1),lattice(indicesinternalC2(1)-100,2),'Marker','pentagram')
        plot(lattice(indicesinternalC3(1)-100,1),lattice(indicesinternalC3(1)-100,2),'Marker','pentagram')
        plot(lattice(indicesinternalC4(1)-100,1),lattice(indicesinternalC4(1)-100,2),'Marker','pentagram') 
        
        plot(lattice(indicesinternalC1(2),1),lattice(indicesinternalC1(2),2),'Marker','square')
        plot(lattice(indicesinternalC2(2),1),lattice(indicesinternalC2(2),2),'Marker','square')
        plot(lattice(indicesinternalC3(2),1),lattice(indicesinternalC3(2),2),'Marker','square')
        plot(lattice(indicesinternalC4(2),1),lattice(indicesinternalC4(2),2),'Marker','square')

        plot(lattice(indicesinternalC1(3),1),lattice(indicesinternalC1(3),2),'Marker','diamond')
        plot(lattice(indicesinternalC2(3),1),lattice(indicesinternalC2(3),2),'Marker','diamond')
        plot(lattice(indicesinternalC3(3),1),lattice(indicesinternalC3(3),2),'Marker','diamond')
        plot(lattice(indicesinternalC4(3),1),lattice(indicesinternalC4(3),2),'Marker','diamond')






%%






flagperiodicity = 0;
periodicity = 0;
flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(logfullfile,N,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);

[metriccoefficients,g,sqrtg] = computemetriccoefficients2D(covariantbase);

contravariantbase = computecontravariantbase2D(covariantbase,sqrtg);

[reciprocalmetriccoefficients,g,sqrtg] = computereciprocalmetriccoefficients2D(contravariantbase);

firstChristoffelsymbol = computefirstChristoffelsymbol2D(N,deltaq,covariantbase,firstdevneighbours);

secondChristoffelsymbol = computesecondChristoffelsymbol2D(N,deltaq,covariantbase,contravariantbase,firstdevneighbours);

Riemanntensor = computeRiemanntensor2D(N,deltaq,secondChristoffelsymbol,metriccoefficients,firstdevneighbours);

Riccitensor = computeRiccitensor2D(Riemanntensor,reciprocalmetriccoefficients);

R = computeRiccicurvature2D(Riccitensor,reciprocalmetriccoefficients);

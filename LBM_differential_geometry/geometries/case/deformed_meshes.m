clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);

flagperiodicity = 0;
periodicity = 0;
flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(N,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);

[metriccoefficients,g,sqrtg] = computemetriccoefficients2D(covariantbase);

contravariantbase = computecontravariantbase2D(covariantbase,sqrtg);

[reciprocalmetriccoefficients,g,sqrtg] = computereciprocalmetriccoefficients2D(contravariantbase);

firstChristoffelsymbol = computefirstChristoffelsymbol2D(N,deltaq,covariantbase,firstdevneighbours);

secondChristoffelsymbol = computesecondChristoffelsymbol2D(N,deltaq,covariantbase,contravariantbase,firstdevneighbours);

Riemanntensor = computeRiemanntensor2D(N,deltaq,secondChristoffelsymbol,metriccoefficients,firstdevneighbours);

Riccitensor = computeRiccitensor2D(Riemanntensor,reciprocalmetriccoefficients);

R = computeRiccicurvature2D(Riccitensor,reciprocalmetriccoefficients);

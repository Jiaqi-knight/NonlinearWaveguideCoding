function[circle]=circularsection2D(x0,y0,r,eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,ymin,ymax,itmax,tol)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 15th, 2014
%    Last update: July 23rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Nx = Neta1;
Ny = Neta2;

eta1max = eta1min + (Neta1-1)*deltaeta1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

flagperiodicity = 0;
periodicity = 0;

flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;

deltaq = [deltaeta1 deltaeta2];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

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

circle = transfiniteinterpolation2D(Nx*Ny,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);

[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(Nx*Ny,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

if itmax~=0
    circle = sparseellipticgridgen2D(Nx,Nx*Ny,circle,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 
end

return
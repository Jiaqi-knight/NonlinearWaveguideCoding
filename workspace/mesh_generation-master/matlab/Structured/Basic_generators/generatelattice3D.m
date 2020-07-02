function[lattice]=generatelattice3D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,eta3min,Neta3,deltaeta3,xmin,xmax,Nx,ymin,ymax,Ny,zmin,zmax,Nz)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Z眉rich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: June 26th, 2014
%    Last update: June 26th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

lattice = zeros(Nx*Ny*Nz,9);%后面为了计算gij等，用到differential geometry

deltax = (xmax - xmin)/(Nx-1);
deltay = (ymax - ymin)/(Ny-1);
deltaz = (zmax - zmin)/(Nz-1);

xs = (xmin:deltax:xmax)';
ys = zeros(Nx*Ny,2);

eta1s = (eta1min:deltaeta1:eta1min+deltaeta1*(Neta1-1))';
eta2s = zeros(Neta1*Neta2,2);

for j=1:Neta2
    eta2s((j-1)*Neta1+1:j*Neta1,1) = eta1s;
    eta2s((j-1)*Neta1+1:j*Neta1,2) = (eta2min+deltaeta2*(j-1))*ones(Neta1,1);
end

clear eta1s

for k=1:Neta3
    lattice((k-1)*Neta1*Neta2+1:k*Neta1*Neta2,1:2) = eta2s;
    lattice((k-1)*Neta1*Neta2+1:k*Neta1*Neta2,3) = (eta3min+deltaeta3*(k-1))*ones(Neta1*Neta2,1);
end

clear eta2s

for j=1:Ny
    ys((j-1)*Nx+1:j*Nx,1) = xs;
    ys((j-1)*Nx+1:j*Nx,2) = (ymin+deltay*(j-1))*ones(Nx,1);
end

clear xs

for k=1:Nz
    lattice((k-1)*Nx*Ny+1:k*Nx*Ny,4:5) = ys;
    lattice((k-1)*Nx*Ny+1:k*Nx*Ny,6) = (zmin+deltaz*(k-1))*ones(Nx*Ny,1);
end

clear ys

lattice(:,7:9) = lattice(:,4:6);

% indeces for finite difference method
% indices = (1:Nx*Ny*Nz)';
% 
% indecesfacexz1 = zeros(Nx*Nz,1);
% indecesfacexz2 = zeros(Nx*Nz,1);
% indecesfacexy1 = zeros(Nx*Ny,1);
% indecesfacexy2 = zeros(Nx*Ny,1);
% indecesfaceyz1 = zeros(Ny*Nz,1);
% indecesfaceyz2 = zeros(Ny*Nz,1);
% 
% for i=1:Nx
%     for k=1:Nz
%         j1 = 1;
%         j2 = Ny;
%         indecesfacexz1(i+(k-1)*Nx) = i+(j1-1)*Nx+(k-1)*Nx*Ny;
%         indecesfacexz2(i+(k-1)*Nx) = i+(j2-1)*Nx+(k-1)*Nx*Ny;
%     end
% end
% for i=1:Nx
%     for j=1:Ny
%         k1 = 1;
%         k2 = Nz;
%         indecesfacexy1(i+(j-1)*Nx) = i+(j-1)*Nx+(k1-1)*Nx*Ny;
%         indecesfacexy2(i+(j-1)*Nx) = i+(j-1)*Nx+(k2-1)*Nx*Ny;
%     end
% end
% for j=1:Ny
%     for k=1:Nz
%         i1 = 1;
%         i2 = Nx;
%         indecesfaceyz1(j+(k-1)*Ny) = i1+(j-1)*Nx+(k-1)*Nx*Ny;
%         indecesfaceyz2(j+(k-1)*Ny) = i2+(j-1)*Nx+(k-1)*Nx*Ny;
%     end
% end
% 
% lattice(:,4) = deltaeta1*ones(Nx*Ny*Nz,1);
% lattice(:,5) = indices - 1;
% lattice(indecesfaceyz1,5) = indecesfaceyz1;
% lattice(:,6) = indices + 1;
% lattice(indecesfaceyz2,6) = 
% lattice(:,7) = deltaeta2*ones(Nx*Ny*Nz,1);
% lattice(:,8) = indices - Nx;
% lattice(:,9) = indices + Nx;
% lattice(:,10) = deltaeta3*ones(Nx*Ny*Nz,1);
% lattice(:,11) = indices - Nx*Ny;
% lattice(:,12) = indices + Nx*Ny;



return
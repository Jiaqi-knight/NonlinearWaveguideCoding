function[lattice]=generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 10th, 2014
%    Last update: July 10th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

lattice = zeros(Nx*Ny,6);

deltax = (xmax - xmin)/(Nx-1);
deltay = (ymax - ymin)/(Ny-1);

xs = (xmin:deltax:xmax)';

eta1s = (eta1min:deltaeta1:eta1min+deltaeta1*(Neta1-1))';

for j=1:Neta2
    lattice((j-1)*Neta1+1:j*Neta1,1) = eta1s;
    lattice((j-1)*Neta1+1:j*Neta1,2) = (eta2min+deltaeta2*(j-1))*ones(Neta1,1);
end

clear eta1s

for j=1:Ny
    lattice((j-1)*Nx+1:j*Nx,3) = xs;
    lattice((j-1)*Nx+1:j*Nx,4) = (ymin+deltay*(j-1))*ones(Nx,1);
end

clear xs

lattice(:,5:6) = lattice(:,3:4);

return
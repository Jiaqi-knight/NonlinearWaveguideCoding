function[cells]=generatecells3D(Nx,Ny,Nz)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 3rd, 2014
%    Last update: July 3rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

cells = zeros((Nx-1)*(Ny-1)*(Nz-1),9);

for i=1:Nz-1
    for j=1:Ny-1
        for k=1:Nx-1
            nodeindex1 = k + (j-1)*Nx + (i-1)*Nx*Ny;
            nodeindex2 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx;
            nodeindex3 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx + 1;
            nodeindex4 = k + (j-1)*Nx + (i-1)*Nx*Ny + 1;
            nodeindex5 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx*Ny;
            nodeindex6 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx*Ny + Nx;
            nodeindex7 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx*Ny + Nx + 1;
            nodeindex8 = k + (j-1)*Nx + (i-1)*Nx*Ny + Nx*Ny + 1;
       
            cellindex = k + (j-1)*(Nx-1) + (i-1)*(Nx-1)*(Ny-1);
                        
            cells(cellindex,1) = 8;
            cells(cellindex,2) = nodeindex1;
            cells(cellindex,3) = nodeindex4;
            cells(cellindex,4) = nodeindex3;
            cells(cellindex,5) = nodeindex2;
            cells(cellindex,6) = nodeindex5;
            cells(cellindex,7) = nodeindex8;
            cells(cellindex,8) = nodeindex7;
            cells(cellindex,9) = nodeindex6;
        end
    end
end

return
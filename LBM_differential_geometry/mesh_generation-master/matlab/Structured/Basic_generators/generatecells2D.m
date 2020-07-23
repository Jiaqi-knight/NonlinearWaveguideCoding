function[cells]=generatecells2D(Nx,Ny)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 30th, 2014
%    Last update: July 30th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

cells = zeros((Nx-1)*(Ny-1),5);

for j=1:Ny-1
    for k=1:Nx-1
        nodeindex1 = k + (j-1)*Nx;
        nodeindex2 = k + (j-1)*Nx + Nx;
        nodeindex3 = k + (j-1)*Nx + Nx + 1;
        nodeindex4 = k + (j-1)*Nx + 1;
       
        cellindex = k + (j-1)*(Nx-1);
                        
        cells(cellindex,1) = 4;
        cells(cellindex,2) = nodeindex1;
        cells(cellindex,3) = nodeindex4;
        cells(cellindex,4) = nodeindex3;
        cells(cellindex,5) = nodeindex2;
    end
end

return
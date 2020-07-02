function[LBMdata]=generate_LBMgeometricinputdata(D,lattice,mesh)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: June 12th, 2014
%    Last update: June 12th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

if lattice.totN==mesh.totN
    switch D
        case 2
            LBMdata = zeros(lattice.totN,4);
            for i=1:lattice.totN
                LBMdata(i,1) = lattice.nodes.meshcoordinates(1,i);
                LBMdata(i,2) = lattice.nodes.meshcoordinates(2,i);
                LBMdata(i,3) = mesh.nodes.meshcoordinates(1,i);
                LBMdata(i,4) = mesh.nodes.meshcoordinates(2,i);
            end
        case 3
            LBMdata = zeros(lattice.totN,6);
            for i=1:lattice.totN
                LBMdata(i,1) = lattice.nodes.meshcoordinates(1,i);
                LBMdata(i,2) = lattice.nodes.meshcoordinates(2,i);
                LBMdata(i,3) = lattice.nodes.meshcoordinates(3,i);
                LBMdata(i,4) = mesh.nodes.meshcoordinates(1,i);
                LBMdata(i,5) = mesh.nodes.meshcoordinates(2,i);
                LBMdata(i,6) = mesh.nodes.meshcoordinates(3,i);
            end
    end
else
    disp('Error: lattice has a different number of nodes from physical mesh. A one-to-one correspondence cannot be established.')
end

return
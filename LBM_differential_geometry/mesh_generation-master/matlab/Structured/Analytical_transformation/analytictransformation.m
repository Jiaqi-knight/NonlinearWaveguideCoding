function[lattice]=analytictransformation(D,lattice,fhandle)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 27th, 2014
%    Last update: August 6th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

switch D
    case 2
        lattice(:,3:4) = fhandle(lattice(:,3),lattice(:,4));
        lattice(:,5:6) = lattice(:,3:4);
    case 3
        lattice(:,4:6) = fhandle(lattice(:,4),lattice(:,5),lattice(:,6));
        lattice(:,7:9) = lattice(:,4:6);  
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return
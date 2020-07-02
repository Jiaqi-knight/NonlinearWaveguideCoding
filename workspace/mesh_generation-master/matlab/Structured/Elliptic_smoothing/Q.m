function[Qvec]=Q(lattice)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 11th, 2014
%    Last update: July 11th, 2014
%
%          Input: meshed computational domain compdomain
%         Output: mesh in the physical domain

%%

Qvec = zeros(size(lattice,1),1);

return
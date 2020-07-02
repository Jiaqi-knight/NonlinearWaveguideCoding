function[edges,triangles]=Delaunaytriangulation2D(p,printflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 23rd, 2014
%    Last update: April 24th, 2014
%
%          Input: N x 4 vector points in space (indices + coordinates)
%         Output: M x 1 vector of indeces of points belonging to the hull

%% Initialize triangulation - Radial sweep

[edges,triangles] = initializeDT2D(p,printflag);

return
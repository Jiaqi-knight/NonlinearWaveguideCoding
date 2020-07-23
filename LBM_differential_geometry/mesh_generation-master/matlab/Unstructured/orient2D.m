function[d]=orient2D(a,b,c)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 15th, 2014
%    Last update: April 15th, 2014
%
%          Input: 2 x 1 vector a of first point coordinates
%                 2 x 1 vector b of second point coordinates
%                 2 x 1 vector c of third point coordinates
%         Output: scalar determinant d, corresponding to the signed area of
%                 the parallelogram determined by the vectors a-c and b-c
%%

d = det([a(1,1)-c(1,1) a(2,1)-c(2,1);...
         b(1,1)-c(1,1) b(2,1)-c(2,1)]);

return
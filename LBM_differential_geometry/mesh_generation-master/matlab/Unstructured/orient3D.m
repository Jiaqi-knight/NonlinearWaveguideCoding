function[e]=orient3D(a,b,c,d)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 15th, 2014
%    Last update: April 16th, 2014
%
%          Input: 3 x 1 vector a of first point coordinates
%                 3 x 1 vector b of second point coordinates
%                 3 x 1 vector c of third point coordinates
%                 3 x 1 vector d of fourth point coordinates
%         Output: scalar determinant e, corresponding to the signed area of
%                 the parallelepiped determined by the vectors a-d, b-d and
%                 c-d                 
%%

e = det([a(1,1)-d(1,1) a(2,1)-d(2,1) a(3,1)-d(3,1);...
         b(1,1)-d(1,1) b(2,1)-d(2,1) b(3,1)-d(3,1);...
         c(1,1)-d(1,1) c(2,1)-d(2,1) c(3,1)-d(3,1)]);

return
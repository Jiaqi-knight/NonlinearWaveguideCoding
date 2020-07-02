function[e]=incircle(a,b,c,d)

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
%          Input: 2 x 1 vector a of first point coordinates
%                 2 x 1 vector b of second point coordinates
%                 2 x 1 vector c of third point coordinates
%                 2 x 1 vector d of fourth point coordinates
%         Output: scalar determinant e, which results to be negative, zero
%                 or positive respectively if point d is inside, on (or 
%                 aligned) or outside the (possibly degenerate) circle 
%                 passing through points a, b and c (in counterclockwise 
%                 order around the circle)
%          Notes: related to the parabolic lifting map
%%

e = sign(det([a(1,1)-d(1,1) a(2,1)-d(2,1) (a(1,1)-d(1,1))^2+(a(2,1)-d(2,1))^2;...
         b(1,1)-d(1,1) b(2,1)-d(2,1) (b(1,1)-d(1,1))^2+(b(2,1)-d(2,1))^2;...
         c(1,1)-d(1,1) c(2,1)-d(2,1) (c(1,1)-d(1,1))^2+(c(2,1)-d(2,1))^2]));

return
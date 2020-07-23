function[f]=insphere(a,b,c,d,e)

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
%          Input: 3 x 1 vector a of first point coordinates
%                 3 x 1 vector b of second point coordinates
%                 3 x 1 vector c of third point coordinates
%                 3 x 1 vector d of fourth point coordinates
%                 3 x 1 vector e of fifth point coordinates
%         Output: scalar determinant f, which results to be negative, zero
%                 or positive respectively if point e is inside, on (or 
%                 coplanar) or outside the (possibly degenerate) sphere 
%                 passing through points a, b, c and d (oriented in such 
%                 a way that orient3D is non-negative)
%          Notes: related to the parabolic lifting map
%%

f = sign(det([a(1,1)-e(1,1) a(2,1)-e(2,1) a(3,1)-e(3,1) (a(1,1)-e(1,1))^2+(a(2,1)-e(2,1))^2+(a(3,1)-e(3,1))^2;...
         b(1,1)-e(1,1) b(2,1)-e(2,1) b(3,1)-e(3,1) (b(1,1)-e(1,1))^2+(b(2,1)-e(2,1))^2+(b(3,1)-e(3,1))^2;...
         b(1,1)-e(1,1) c(2,1)-e(2,1) c(3,1)-e(3,1) (c(1,1)-e(1,1))^2+(c(2,1)-e(2,1))^2+(c(3,1)-e(3,1))^2;...
         d(1,1)-e(1,1) d(2,1)-e(2,1) d(3,1)-e(3,1) (d(1,1)-e(1,1))^2+(d(2,1)-e(2,1))^2+(d(3,1)-e(3,1))^2]));

return
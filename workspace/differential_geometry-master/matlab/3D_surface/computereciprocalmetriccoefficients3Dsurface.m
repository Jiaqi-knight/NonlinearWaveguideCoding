function[reciprocalmetriccoefficients,g,sqrtg]=computereciprocalmetriccoefficients3Dsurface(contravariantbase)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 17th, 2014
%    Last update: July 17th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

g11 = sum(contravariantbase(:,1:3).*contravariantbase(:,1:3),2);

g22 = sum(contravariantbase(:,4:6).*contravariantbase(:,4:6),2);

g12 = sum(contravariantbase(:,1:3).*contravariantbase(:,4:6),2);

g = g11.*g22 - g12.*g12;

sqrtg = sqrt(g);

reciprocalmetriccoefficients = [g11 g22 g12];

return
function[reciprocalmetriccoefficients,g,sqrtg]=computereciprocalmetriccoefficients3D(contravariantbase)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 28th, 2014
%    Last update: July 8th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

g11 = sum(contravariantbase(:,1:3).*contravariantbase(:,1:3),2);

g22 = sum(contravariantbase(:,4:6).*contravariantbase(:,4:6),2);

g33 = sum(contravariantbase(:,7:9).*contravariantbase(:,7:9),2);

g12 = sum(contravariantbase(:,1:3).*contravariantbase(:,4:6),2);

g13 = sum(contravariantbase(:,1:3).*contravariantbase(:,7:9),2);

g23 = sum(contravariantbase(:,4:6).*contravariantbase(:,7:9),2);

g = g11.*(g22.*g33-g23.^2) + g12.*(g13.*g23-g12.*g33) + g13.*(g12.*g23-g13.*g22);

sqrtg = sqrt(g);

reciprocalmetriccoefficients = [g11 g22 g33 g12 g13 g23];

return
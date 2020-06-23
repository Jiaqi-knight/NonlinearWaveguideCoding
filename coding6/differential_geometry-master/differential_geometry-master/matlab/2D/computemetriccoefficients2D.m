function[metriccoefficients,g,sqrtg]=computemetriccoefficients2D(covariantbase)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 11th, 2014
%    Last update: July 11th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

g11 = sum(covariantbase(:,1:2).*covariantbase(:,1:2),2);

g22 = sum(covariantbase(:,3:4).*covariantbase(:,3:4),2);

g12 = sum(covariantbase(:,1:2).*covariantbase(:,3:4),2);

g = g11.*g22 - g12.*g12;

sqrtg = sqrt(g);

metriccoefficients = [g11 g22 g12];

return
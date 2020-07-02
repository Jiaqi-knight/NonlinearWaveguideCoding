function[coords]=sphericalcoordinates(r,phi,theta)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: August 6th, 2014
%    Last update: August 6th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

coords = [r.*sin(phi).*cos(theta) r.*sin(phi).*sin(theta) r.*cos(phi)];

return
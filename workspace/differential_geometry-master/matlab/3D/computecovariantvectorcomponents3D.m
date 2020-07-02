function[vi]=computecovariantvectorcomponents3D(covariantbase,v)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 28th, 2014
%    Last update: May 28th, 2014
%
%    Description: 
%          Input: N x 3 vector v
%         Output: 

%%

g1 = covariantbase(:,1:3);
g2 = covariantbase(:,4:6);
g3 = covariantbase(:,7:9);

v1 = sum(v.*g1,2);
v2 = sum(v.*g2,2);
v3 = sum(v.*g3,2);

vi = [v1 v2 v3];

return
function[vI]=computecontravariantvectorcomponents3D(contravariantbase)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 28th, 2014
%    Last update: May 30th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

G1 = contravariantbase(:,1:3);
G2 = contravariantbase(:,4:6);
G3 = contravariantbase(:,7:9);

v1 = sum(v.*G1,2);
v2 = sum(v.*G2,2);
v3 = sum(v.*G3,2);

vI = [v1 v2 v3];

return
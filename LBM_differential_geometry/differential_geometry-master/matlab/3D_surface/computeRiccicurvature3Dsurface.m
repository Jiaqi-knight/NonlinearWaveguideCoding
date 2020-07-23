function[R]=computeRiccicurvature3Dsurface(Riccitensor,reciprocalmetriccoefficients)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 18th, 2014
%    Last update: July 18th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

G11 = reciprocalmetriccoefficients(:,1);
G22 = reciprocalmetriccoefficients(:,2);
G12 = reciprocalmetriccoefficients(:,3);

R11 = Riccitensor(:,1);
R22 = Riccitensor(:,2);
R12 = Riccitensor(:,3);

R = G11*R11 + G22*R22 + G12*R12;

return
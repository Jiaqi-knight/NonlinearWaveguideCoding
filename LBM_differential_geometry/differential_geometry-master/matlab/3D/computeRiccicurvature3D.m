function[R]=computeRiccicurvature3D(Riccitensor,reciprocalmetriccoefficients)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 29th, 2014
%    Last update: May 29th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

G11 = reciprocalmetriccoefficients(:,1);
G22 = reciprocalmetriccoefficients(:,2);
G33 = reciprocalmetriccoefficients(:,3);
G12 = reciprocalmetriccoefficients(:,4);
G13 = reciprocalmetriccoefficients(:,5);
G23 = reciprocalmetriccoefficients(:,6);

R11 = Riccitensor(:,1);
R22 = Riccitensor(:,2);
R33 = Riccitensor(:,3);
R12 = Riccitensor(:,4);
R13 = Riccitensor(:,5);
R23 = Riccitensor(:,6);

R = G11.*R11 + G22.*R22 + G33.*R33 + G12.*R12 + G13.*R13 + G23.*R23;

return
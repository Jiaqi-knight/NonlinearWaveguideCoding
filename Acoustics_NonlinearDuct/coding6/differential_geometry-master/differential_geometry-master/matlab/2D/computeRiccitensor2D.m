function[Riccitensor]=computeRiccitensor2D(Riemanntensor,reciprocalmetriccoefficients)

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

R1111 = Riemanntensor(:,1);
R1112 = Riemanntensor(:,2);
R1121 = Riemanntensor(:,3);
R1122 = Riemanntensor(:,4);
R1211 = Riemanntensor(:,5);
R1212 = Riemanntensor(:,6);
R1221 = Riemanntensor(:,7);
R1222 = Riemanntensor(:,8);
%R2111 = Riemanntensor(:,9);
%R2112 = Riemanntensor(:,10);
R2121 = Riemanntensor(:,11);
R2122 = Riemanntensor(:,12);
%R2211 = Riemanntensor(:,13);
%R2212 = Riemanntensor(:,14);
R2221 = Riemanntensor(:,15);
R2222 = Riemanntensor(:,16);

% R_ij = R_ji
R11 = G11.*R1111 + G12.*R1112 + G12.*R1211 + G22.*R1212;
R12 = G11.*R1121 + G12.*R1122 + G12.*R1221 + G22.*R1222;
R22 = G11.*R2121 + G12.*R2122 + G12.*R2221 + G22.*R2222;

Riccitensor = [R11 R22 R12];

return
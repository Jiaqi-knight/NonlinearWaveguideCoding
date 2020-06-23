function[contravariantbase]=computecontravariantbase3D(covariantbase,sqrtg)

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
%          Input: 
%         Output: 

%%

g1 = covariantbase(:,1:3);
g2 = covariantbase(:,4:6);
g3 = covariantbase(:,7:9);

G1 = [(g2(:,2).*g3(:,3)-g2(:,3).*g3(:,2)) (g2(:,3).*g3(:,1)-g2(:,1).*g3(:,3)) (g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1))]./[sqrtg sqrtg sqrtg];

G2 = [(g3(:,2).*g1(:,3)-g3(:,3).*g1(:,2)) (g3(:,3).*g1(:,1)-g3(:,1).*g1(:,3)) (g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1))]./[sqrtg sqrtg sqrtg];

G3 = [(g1(:,2).*g2(:,3)-g1(:,3).*g2(:,2)) (g1(:,3).*g2(:,1)-g1(:,1).*g2(:,3)) (g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1))]./[sqrtg sqrtg sqrtg];

contravariantbase = [G1 G2 G3];

return
function[contravariantbase]=computecontravariantbase2D(covariantbase,sqrtg)

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

N = size(covariantbase,1);

g1 = [covariantbase(:,1:2) zeros(N,1)];
g2 = [covariantbase(:,3:4) zeros(N,1)];

vecprod = [zeros(N,1) zeros(N,1) (g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1))];
norm = sqrt(sum(vecprod.^2,2));
g3 = vecprod./[norm norm norm];

G1 = [g2(:,2).*g3(:,3)-g2(:,3).*g3(:,2) g2(:,3).*g3(:,1)-g2(:,1).*g3(:,3)]./[sqrtg sqrtg]; %g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1)]./sqrtg;

G2 = [g3(:,2).*g1(:,3)-g3(:,3).*g1(:,2) g3(:,3).*g1(:,1)-g3(:,1).*g1(:,3)]./[sqrtg sqrtg]; %g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1)]./sqrtg;

%G3 = [g1(:,2).*g2(:,3)-g1(:,3).*g2(:,2) g1(:,3).*g2(:,1)-g1(:,1).*g2(:,3) g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1)]./sqrtg;

contravariantbase = [G1 G2];

return
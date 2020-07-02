function[covariantbase]=computecovariantbase3Dsurface(N,deltaq,lattice,firstdevneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 17th, 2014
%    Last update: July 17th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

covariantbase = zeros(N,9);

for i=1:N
    for j=1:2
        switch firstdevneighbours(i,4*(j-1)+1)
            case 1
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = 0.5.*[(lattice(firstdevneighbours(i,4*(j-1)+3),7)-lattice(firstdevneighbours(i,4*(j-1)+2),7)) (lattice(firstdevneighbours(i,4*(j-1)+3),8)-lattice(firstdevneighbours(i,4*(j-1)+2),8)) (lattice(firstdevneighbours(i,4*(j-1)+3),9)-lattice(firstdevneighbours(i,4*(j-1)+2),9))]./deltaq(j);
            case 2
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = [(-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),7)+2*lattice(firstdevneighbours(i,4*(j-1)+3),7)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),7)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),8)+2*lattice(firstdevneighbours(i,4*(j-1)+3),8)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),8)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),9)+2*lattice(firstdevneighbours(i,4*(j-1)+3),9)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),9))]/deltaq(j);
            case 3
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = [(1.5*lattice(firstdevneighbours(i,4*(j-1)+2),7)-2*lattice(firstdevneighbours(i,4*(j-1)+3),7)+0.5*lattice(firstdevneighbours(i,4*(j-1)+4),7)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),8)+2*lattice(firstdevneighbours(i,4*(j-1)+3),8)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),8)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),9)+2*lattice(firstdevneighbours(i,4*(j-1)+3),9)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),9))]/deltaq(j);
        end
    end
    vecprod = [(covariantbase(i,2)*covariantbase(i,6)-covariantbase(i,3)*covariantbase(i,5)) (covariantbase(i,3)*covariantbase(i,4)-covariantbase(i,1)*covariantbase(i,6)) (covariantbase(i,1)*covariantbase(i,5)-covariantbase(i,2)*covariantbase(i,4))];
    norm = sqrt(sum(vecprod.^2,2));
    covariantbase(i,7:9) = vecprod./norm;
end

return
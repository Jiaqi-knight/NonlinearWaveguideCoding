function[covariantbase]=computecovariantbase2D(N,deltaq,lattice,firstdevneighbours)

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

covariantbase = zeros(N,4);

for i=1:N
    for j=1:2
        switch firstdevneighbours(i,4*(j-1)+1)
            case 1
                covariantbase(i,2*(j-1)+1:2*(j-1)+2) = 0.5.*[(lattice(firstdevneighbours(i,4*(j-1)+3),5)-lattice(firstdevneighbours(i,4*(j-1)+2),5)) (lattice(firstdevneighbours(i,4*(j-1)+3),6)-lattice(firstdevneighbours(i,4*(j-1)+2),6))]./deltaq(j);
            case 2
                covariantbase(i,2*(j-1)+1:2*(j-1)+2) = [(-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),5)+2*lattice(firstdevneighbours(i,4*(j-1)+3),5)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),5)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),6)+2*lattice(firstdevneighbours(i,4*(j-1)+3),6)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),6))]/deltaq(j);
            case 3
                covariantbase(i,2*(j-1)+1:2*(j-1)+2) = [(1.5*lattice(firstdevneighbours(i,4*(j-1)+2),5)-2*lattice(firstdevneighbours(i,4*(j-1)+3),5)+0.5*lattice(firstdevneighbours(i,4*(j-1)+4),5)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),6)+2*lattice(firstdevneighbours(i,4*(j-1)+3),6)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),6))]/deltaq(j);
        end
    end
end
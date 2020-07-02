function[covariantbase]=computecovariantbase3D(N,deltaq,lattice,firstdevneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 28th, 2014
%    Last update: July 8th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

covariantbase = zeros(N,9);

for i=1:N
    for j=1:3
        switch firstdevneighbours(i,4*(j-1)+1)
            case 1
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = 0.5.*[(lattice(firstdevneighbours(i,4*(j-1)+3),7)-lattice(firstdevneighbours(i,4*(j-1)+2),7)) (lattice(firstdevneighbours(i,4*(j-1)+3),8)-lattice(firstdevneighbours(i,4*(j-1)+2),8)) (lattice(firstdevneighbours(i,4*(j-1)+3),9)-lattice(firstdevneighbours(i,4*(j-1)+2),9))]./deltaq(j);
            case 2
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = [(-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),7)+2*lattice(firstdevneighbours(i,4*(j-1)+3),7)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),7)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),8)+2*lattice(firstdevneighbours(i,4*(j-1)+3),8)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),8)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),9)+2*lattice(firstdevneighbours(i,4*(j-1)+3),9)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),9))]/deltaq(j);
            case 3
                covariantbase(i,3*(j-1)+1:3*(j-1)+3) = [(1.5*lattice(firstdevneighbours(i,4*(j-1)+2),7)-2*lattice(firstdevneighbours(i,4*(j-1)+3),7)+0.5*lattice(firstdevneighbours(i,4*(j-1)+4),7)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),8)+2*lattice(firstdevneighbours(i,4*(j-1)+3),8)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),8)) (-1.5*lattice(firstdevneighbours(i,4*(j-1)+2),9)+2*lattice(firstdevneighbours(i,4*(j-1)+3),9)-0.5*lattice(firstdevneighbours(i,4*(j-1)+4),9))]/deltaq(j);
        end
    end
end

% covariantbase = [(geometricdata(:,8).*geometricdata(geometricdata(:,5),1)+geometricdata(:,9).*geometricdata(geometricdata(:,6),1)+geometricdata(:,10).*geometricdata(geometricdata(:,7),1))./geometricdata(:,4)...
%                  (geometricdata(:,8).*geometricdata(geometricdata(:,5),2)+geometricdata(:,9).*geometricdata(geometricdata(:,6),2)+geometricdata(:,10).*geometricdata(geometricdata(:,7),2))./geometricdata(:,4)...
%                  (geometricdata(:,8).*geometricdata(geometricdata(:,5),3)+geometricdata(:,9).*geometricdata(geometricdata(:,6),3)+geometricdata(:,10).*geometricdata(geometricdata(:,7),3))./geometricdata(:,4)...
%                  (geometricdata(:,15).*geometricdata(geometricdata(:,12),1)+geometricdata(:,16).*geometricdata(geometricdata(:,13),1)+geometricdata(:,17).*geometricdata(geometricdata(:,14),1))./geometricdata(:,11)...
%                  (geometricdata(:,15).*geometricdata(geometricdata(:,12),2)+geometricdata(:,16).*geometricdata(geometricdata(:,13),2)+geometricdata(:,17).*geometricdata(geometricdata(:,14),2))./geometricdata(:,11)...
%                  (geometricdata(:,15).*geometricdata(geometricdata(:,12),3)+geometricdata(:,16).*geometricdata(geometricdata(:,13),3)+geometricdata(:,17).*geometricdata(geometricdata(:,14),3))./geometricdata(:,11)...
%                  (geometricdata(:,22).*geometricdata(geometricdata(:,19),1)+geometricdata(:,23).*geometricdata(geometricdata(:,20),1)+geometricdata(:,24).*geometricdata(geometricdata(:,21),1))./geometricdata(:,18)...
%                  (geometricdata(:,22).*geometricdata(geometricdata(:,19),2)+geometricdata(:,23).*geometricdata(geometricdata(:,20),2)+geometricdata(:,24).*geometricdata(geometricdata(:,21),2))./geometricdata(:,18)...
%                  (geometricdata(:,22).*geometricdata(geometricdata(:,19),3)+geometricdata(:,23).*geometricdata(geometricdata(:,20),3)+geometricdata(:,24).*geometricdata(geometricdata(:,21),3))./geometricdata(:,18)];

return
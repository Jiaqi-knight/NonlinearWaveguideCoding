function[firstChristoffelsymbol]=computefirstChristoffelsymbol2D(N,deltaq,covariantbase,firstdevneighbours)

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

g1 = covariantbase(:,1:2);
g2 = covariantbase(:,3:4);
%g3 = covariantbase(:,7:9);

g1d1 = zeros(N,2);
g1d2 = zeros(N,2);
g2d1 = zeros(N,2);
g2d2 = zeros(N,2);

for i=1:N
    j1 = 1;
    switch firstdevneighbours(i,4*(j1-1)+1)
        case 1
            g1d1(i,1:2) = 0.5.*[(g1(firstdevneighbours(i,4*(j1-1)+3),1)-g1(firstdevneighbours(i,4*(j1-1)+2),1)) (g1(firstdevneighbours(i,4*(j1-1)+3),2)-g1(firstdevneighbours(i,4*(j1-1)+2),2))]./deltaq(j1);
            g2d1(i,1:2) = 0.5.*[(g2(firstdevneighbours(i,4*(j1-1)+3),1)-g2(firstdevneighbours(i,4*(j1-1)+2),1)) (g2(firstdevneighbours(i,4*(j1-1)+3),2)-g2(firstdevneighbours(i,4*(j1-1)+2),2))]./deltaq(j1);
        case 2
            g1d1(i,1:2) = [(-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),1)+2*g1(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),2)+2*g1(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),2))]./deltaq(j1);
            g2d1(i,1:2) = [(-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),1)+2*g2(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),2)+2*g2(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),2))]./deltaq(j1);
        case 3
            g1d1(i,1:2) = [(1.5*g1(firstdevneighbours(i,4*(j1-1)+2),1)-2*g1(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*g1(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),2)+2*g1(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),2))]./deltaq(j1);
            g2d1(i,1:2) = [(1.5*g2(firstdevneighbours(i,4*(j1-1)+2),1)-2*g2(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*g2(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),2)+2*g2(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),2))]./deltaq(j1);
    end
    j2 = 2;
    switch firstdevneighbours(i,4*(j2-1)+1)
        case 1
            g1d2(i,1:2) = 0.5.*[(g1(firstdevneighbours(i,4*(j2-1)+3),1)-g1(firstdevneighbours(i,4*(j2-1)+2),1)) (g1(firstdevneighbours(i,4*(j2-1)+3),2)-g1(firstdevneighbours(i,4*(j2-1)+2),2))]./deltaq(j2);
            g2d2(i,1:2) = 0.5.*[(g2(firstdevneighbours(i,4*(j2-1)+3),1)-g2(firstdevneighbours(i,4*(j2-1)+2),1)) (g2(firstdevneighbours(i,4*(j2-1)+3),2)-g2(firstdevneighbours(i,4*(j2-1)+2),2))]./deltaq(j2);
        case 2
            g1d2(i,1:2) = [(-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),1)+2*g1(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),2)+2*g1(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),2))]./deltaq(j2);
            g2d2(i,1:2) = [(-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),1)+2*g2(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),2)+2*g2(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),2))]./deltaq(j2);
        case 3
            g1d2(i,1:2) = [(1.5*g1(firstdevneighbours(i,4*(j2-1)+2),1)-2*g1(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*g1(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),2)+2*g1(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),2))]./deltaq(j2);
            g2d2(i,1:2) = [(1.5*g2(firstdevneighbours(i,4*(j2-1)+2),1)-2*g2(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*g2(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),2)+2*g2(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),2))]./deltaq(j2);
            
    end
end

firstChristoffelsymbol = [sum(g1d1.*g1,2) sum(g1d2.*g1,2) sum(g2d1.*g1,2) sum(g2d2.*g1,2) ...
                          sum(g1d1.*g2,2) sum(g1d2.*g2,2) sum(g2d1.*g2,2) sum(g2d2.*g2,2) ];

return

function[secondChristoffelsymbol]=computesecondChristoffelsymbol3D(N,deltaq,covariantbase,contravariantbase,firstdevneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 29th, 2014
%    Last update: July 8th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

g1 = covariantbase(:,1:3);
g2 = covariantbase(:,4:6);
g3 = covariantbase(:,7:9);

G1 = contravariantbase(:,1:3);
G2 = contravariantbase(:,4:6);
G3 = contravariantbase(:,7:9);

g1d1 = zeros(N,3);
g1d2 = zeros(N,3);
g1d3 = zeros(N,3);
g2d1 = zeros(N,3);
g2d2 = zeros(N,3);
g2d3 = zeros(N,3);
g3d1 = zeros(N,3);
g3d2 = zeros(N,3);
g3d3 = zeros(N,3);

for i=1:N
    j1 = 1;
    switch firstdevneighbours(i,4*(j1-1)+1)
        case 1
            g1d1(i,1:3) = 0.5.*[(g1(firstdevneighbours(i,4*(j1-1)+3),1)-g1(firstdevneighbours(i,4*(j1-1)+2),1)) (g1(firstdevneighbours(i,4*(j1-1)+3),2)-g1(firstdevneighbours(i,4*(j1-1)+2),2)) (g1(firstdevneighbours(i,4*(j1-1)+3),3)-g1(firstdevneighbours(i,4*(j1-1)+2),3))]./deltaq(j1);
            g2d1(i,1:3) = 0.5.*[(g2(firstdevneighbours(i,4*(j1-1)+3),1)-g2(firstdevneighbours(i,4*(j1-1)+2),1)) (g2(firstdevneighbours(i,4*(j1-1)+3),2)-g2(firstdevneighbours(i,4*(j1-1)+2),2)) (g2(firstdevneighbours(i,4*(j1-1)+3),3)-g2(firstdevneighbours(i,4*(j1-1)+2),3))]./deltaq(j1);
            g3d1(i,1:3) = 0.5.*[(g3(firstdevneighbours(i,4*(j1-1)+3),1)-g3(firstdevneighbours(i,4*(j1-1)+2),1)) (g3(firstdevneighbours(i,4*(j1-1)+3),2)-g3(firstdevneighbours(i,4*(j1-1)+2),2)) (g3(firstdevneighbours(i,4*(j1-1)+3),3)-g3(firstdevneighbours(i,4*(j1-1)+2),3))]./deltaq(j1);
        case 2
            g1d1(i,1:3) = [(-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),1)+2*g1(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),2)+2*g1(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),3)+2*g1(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
            g2d1(i,1:3) = [(-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),1)+2*g2(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),2)+2*g2(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),3)+2*g2(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
            g3d1(i,1:3) = [(-1.5*g3(firstdevneighbours(i,4*(j1-1)+2),1)+2*g3(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*g3(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j1-1)+2),2)+2*g3(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j1-1)+2),3)+2*g3(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
        case 3
            g1d1(i,1:3) = [(1.5*g1(firstdevneighbours(i,4*(j1-1)+2),1)-2*g1(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*g1(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),2)+2*g1(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j1-1)+2),3)+2*g1(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
            g2d1(i,1:3) = [(1.5*g2(firstdevneighbours(i,4*(j1-1)+2),1)-2*g2(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*g2(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),2)+2*g2(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j1-1)+2),3)+2*g2(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
            g3d1(i,1:3) = [(1.5*g3(firstdevneighbours(i,4*(j1-1)+2),1)-2*g3(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*g3(firstdevneighbours(i,4*(j1-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j1-1)+2),2)+2*g3(firstdevneighbours(i,4*(j1-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j1-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j1-1)+2),3)+2*g3(firstdevneighbours(i,4*(j1-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j1-1)+4),3))]/deltaq(j1);
    end
    j2 = 2;
    switch firstdevneighbours(i,4*(j2-1)+1)
        case 1
            g1d2(i,1:3) = 0.5.*[(g1(firstdevneighbours(i,4*(j2-1)+3),1)-g1(firstdevneighbours(i,4*(j2-1)+2),1)) (g1(firstdevneighbours(i,4*(j2-1)+3),2)-g1(firstdevneighbours(i,4*(j2-1)+2),2)) (g1(firstdevneighbours(i,4*(j2-1)+3),3)-g1(firstdevneighbours(i,4*(j2-1)+2),3))]./deltaq(j2);
            g2d2(i,1:3) = 0.5.*[(g2(firstdevneighbours(i,4*(j2-1)+3),1)-g2(firstdevneighbours(i,4*(j2-1)+2),1)) (g2(firstdevneighbours(i,4*(j2-1)+3),2)-g2(firstdevneighbours(i,4*(j2-1)+2),2)) (g2(firstdevneighbours(i,4*(j2-1)+3),3)-g2(firstdevneighbours(i,4*(j2-1)+2),3))]./deltaq(j2);
            g3d2(i,1:3) = 0.5.*[(g3(firstdevneighbours(i,4*(j2-1)+3),1)-g3(firstdevneighbours(i,4*(j2-1)+2),1)) (g3(firstdevneighbours(i,4*(j2-1)+3),2)-g3(firstdevneighbours(i,4*(j2-1)+2),2)) (g3(firstdevneighbours(i,4*(j2-1)+3),3)-g3(firstdevneighbours(i,4*(j2-1)+2),3))]./deltaq(j2);
        case 2
            g1d2(i,1:3) = [(-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),1)+2*g1(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),2)+2*g1(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),3)+2*g1(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
            g2d2(i,1:3) = [(-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),1)+2*g2(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),2)+2*g2(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),3)+2*g2(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
            g3d2(i,1:3) = [(-1.5*g3(firstdevneighbours(i,4*(j2-1)+2),1)+2*g3(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*g3(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j2-1)+2),2)+2*g3(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j2-1)+2),3)+2*g3(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
        case 3
            g1d2(i,1:3) = [(1.5*g1(firstdevneighbours(i,4*(j2-1)+2),1)-2*g1(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*g1(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),2)+2*g1(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j2-1)+2),3)+2*g1(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
            g2d2(i,1:3) = [(1.5*g2(firstdevneighbours(i,4*(j2-1)+2),1)-2*g2(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*g2(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),2)+2*g2(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j2-1)+2),3)+2*g2(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
            g3d2(i,1:3) = [(1.5*g3(firstdevneighbours(i,4*(j2-1)+2),1)-2*g3(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*g3(firstdevneighbours(i,4*(j2-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j2-1)+2),2)+2*g3(firstdevneighbours(i,4*(j2-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j2-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j2-1)+2),3)+2*g3(firstdevneighbours(i,4*(j2-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j2-1)+4),3))]/deltaq(j2);
    end
    j3 = 3;
    switch firstdevneighbours(i,4*(j3-1)+1)
        case 1
            g1d3(i,1:3) = 0.5.*[(g1(firstdevneighbours(i,4*(j3-1)+3),1)-g1(firstdevneighbours(i,4*(j3-1)+2),1)) (g1(firstdevneighbours(i,4*(j3-1)+3),2)-g1(firstdevneighbours(i,4*(j3-1)+2),2)) (g1(firstdevneighbours(i,4*(j3-1)+3),3)-g1(firstdevneighbours(i,4*(j3-1)+2),3))]./deltaq(j3);
            g2d3(i,1:3) = 0.5.*[(g2(firstdevneighbours(i,4*(j3-1)+3),1)-g2(firstdevneighbours(i,4*(j3-1)+2),1)) (g2(firstdevneighbours(i,4*(j3-1)+3),2)-g2(firstdevneighbours(i,4*(j3-1)+2),2)) (g2(firstdevneighbours(i,4*(j3-1)+3),3)-g2(firstdevneighbours(i,4*(j3-1)+2),3))]./deltaq(j3);
            g3d3(i,1:3) = 0.5.*[(g3(firstdevneighbours(i,4*(j3-1)+3),1)-g3(firstdevneighbours(i,4*(j3-1)+2),1)) (g3(firstdevneighbours(i,4*(j3-1)+3),2)-g3(firstdevneighbours(i,4*(j3-1)+2),2)) (g3(firstdevneighbours(i,4*(j3-1)+3),3)-g3(firstdevneighbours(i,4*(j3-1)+2),3))]./deltaq(j3);
        case 2
            g1d3(i,1:3) = [(-1.5*g1(firstdevneighbours(i,4*(j3-1)+2),1)+2*g1(firstdevneighbours(i,4*(j3-1)+3),1)-0.5*g1(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j3-1)+2),2)+2*g1(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j3-1)+2),3)+2*g1(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
            g2d3(i,1:3) = [(-1.5*g2(firstdevneighbours(i,4*(j3-1)+2),1)+2*g2(firstdevneighbours(i,4*(j3-1)+3),1)-0.5*g2(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j3-1)+2),2)+2*g2(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j3-1)+2),3)+2*g2(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
            g3d3(i,1:3) = [(-1.5*g3(firstdevneighbours(i,4*(j3-1)+2),1)+2*g3(firstdevneighbours(i,4*(j3-1)+3),1)-0.5*g3(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j3-1)+2),2)+2*g3(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j3-1)+2),3)+2*g3(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
        case 3
            g1d3(i,1:3) = [(1.5*g1(firstdevneighbours(i,4*(j3-1)+2),1)-2*g1(firstdevneighbours(i,4*(j3-1)+3),1)+0.5*g1(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g1(firstdevneighbours(i,4*(j3-1)+2),2)+2*g1(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g1(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g1(firstdevneighbours(i,4*(j3-1)+2),3)+2*g1(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g1(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
            g2d3(i,1:3) = [(1.5*g2(firstdevneighbours(i,4*(j3-1)+2),1)-2*g2(firstdevneighbours(i,4*(j3-1)+3),1)+0.5*g2(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g2(firstdevneighbours(i,4*(j3-1)+2),2)+2*g2(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g2(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g2(firstdevneighbours(i,4*(j3-1)+2),3)+2*g2(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g2(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
            g3d3(i,1:3) = [(1.5*g3(firstdevneighbours(i,4*(j3-1)+2),1)-2*g3(firstdevneighbours(i,4*(j3-1)+3),1)+0.5*g3(firstdevneighbours(i,4*(j3-1)+4),1)) (-1.5*g3(firstdevneighbours(i,4*(j3-1)+2),2)+2*g3(firstdevneighbours(i,4*(j3-1)+3),2)-0.5*g3(firstdevneighbours(i,4*(j3-1)+4),2)) (-1.5*g3(firstdevneighbours(i,4*(j3-1)+2),3)+2*g3(firstdevneighbours(i,4*(j3-1)+3),3)-0.5*g3(firstdevneighbours(i,4*(j3-1)+4),3))]/deltaq(j3);
    end
end

secondChristoffelsymbol = [sum(g1d1.*G1,2) sum(g1d2.*G1,2) sum(g1d3.*G1,2) sum(g2d1.*G1,2) sum(g2d2.*G1,2) sum(g2d3.*G1,2) sum(g3d1.*G1,2) sum(g3d2.*G1,2) sum(g3d3.*G1,2)...
                           sum(g1d1.*G2,2) sum(g1d2.*G2,2) sum(g1d3.*G2,2) sum(g2d1.*G2,2) sum(g2d2.*G2,2) sum(g2d3.*G2,2) sum(g3d1.*G2,2) sum(g3d2.*G2,2) sum(g3d3.*G2,2)...
                           sum(g1d1.*G3,2) sum(g1d2.*G3,2) sum(g1d3.*G3,2) sum(g2d1.*G3,2) sum(g2d2.*G3,2) sum(g2d3.*G3,2) sum(g3d1.*G3,2) sum(g3d2.*G3,2) sum(g3d3.*G3,2)];

return

% g1d1 = [(geometricdata(:,8).*g1(geometricdata(:,5),1)+geometricdata(:,9).*g1(geometricdata(:,6),1)+geometricdata(:,10).*g1(geometricdata(:,7),1))./geometricdata(:,4)...
%         (geometricdata(:,8).*g1(geometricdata(:,5),2)+geometricdata(:,9).*g1(geometricdata(:,6),2)+geometricdata(:,10).*g1(geometricdata(:,7),2))./geometricdata(:,4)...
%         (geometricdata(:,8).*g1(geometricdata(:,5),3)+geometricdata(:,9).*g1(geometricdata(:,6),3)+geometricdata(:,10).*g1(geometricdata(:,7),3))./geometricdata(:,4)];
%     
% g1d2 = [(geometricdata(:,15).*g1(geometricdata(:,12),1)+geometricdata(:,16).*g1(geometricdata(:,13),1)+geometricdata(:,17).*g1(geometricdata(:,14),1))./geometricdata(:,11)...
%         (geometricdata(:,15).*g1(geometricdata(:,12),2)+geometricdata(:,16).*g1(geometricdata(:,13),2)+geometricdata(:,17).*g1(geometricdata(:,14),2))./geometricdata(:,11)...
%         (geometricdata(:,15).*g1(geometricdata(:,12),3)+geometricdata(:,16).*g1(geometricdata(:,13),3)+geometricdata(:,17).*g1(geometricdata(:,14),3))./geometricdata(:,11)];
% 
% g1d3 = [(geometricdata(:,22).*g1(geometricdata(:,19),1)+geometricdata(:,23).*g1(geometricdata(:,20),1)+geometricdata(:,24).*g1(geometricdata(:,21),1))./geometricdata(:,18)...
%         (geometricdata(:,22).*g1(geometricdata(:,19),2)+geometricdata(:,23).*g1(geometricdata(:,20),2)+geometricdata(:,24).*g1(geometricdata(:,21),2))./geometricdata(:,18)...
%         (geometricdata(:,22).*g1(geometricdata(:,19),3)+geometricdata(:,23).*g1(geometricdata(:,20),3)+geometricdata(:,24).*g1(geometricdata(:,21),3))./geometricdata(:,18)];
%     
% g2d1 = [(geometricdata(:,8).*g2(geometricdata(:,5),1)+geometricdata(:,9).*g2(geometricdata(:,6),1)+geometricdata(:,10).*g2(geometricdata(:,7),1))./geometricdata(:,4)...
%         (geometricdata(:,8).*g2(geometricdata(:,5),2)+geometricdata(:,9).*g2(geometricdata(:,6),2)+geometricdata(:,10).*g2(geometricdata(:,7),2))./geometricdata(:,4)...
%         (geometricdata(:,8).*g2(geometricdata(:,5),3)+geometricdata(:,9).*g2(geometricdata(:,6),3)+geometricdata(:,10).*g2(geometricdata(:,7),3))./geometricdata(:,4)];
%     
% g2d2 = [(geometricdata(:,15).*g2(geometricdata(:,12),1)+geometricdata(:,16).*g2(geometricdata(:,13),1)+geometricdata(:,17).*g2(geometricdata(:,14),1))./geometricdata(:,11)...
%         (geometricdata(:,15).*g2(geometricdata(:,12),2)+geometricdata(:,16).*g2(geometricdata(:,13),2)+geometricdata(:,17).*g2(geometricdata(:,14),2))./geometricdata(:,11)...
%         (geometricdata(:,15).*g2(geometricdata(:,12),3)+geometricdata(:,16).*g2(geometricdata(:,13),3)+geometricdata(:,17).*g2(geometricdata(:,14),3))./geometricdata(:,11)];
% 
% g2d3 = [(geometricdata(:,22).*g2(geometricdata(:,19),1)+geometricdata(:,23).*g2(geometricdata(:,20),1)+geometricdata(:,24).*g2(geometricdata(:,21),1))./geometricdata(:,18)...
%         (geometricdata(:,22).*g2(geometricdata(:,19),2)+geometricdata(:,23).*g2(geometricdata(:,20),2)+geometricdata(:,24).*g2(geometricdata(:,21),2))./geometricdata(:,18)...
%         (geometricdata(:,22).*g2(geometricdata(:,19),3)+geometricdata(:,23).*g2(geometricdata(:,20),3)+geometricdata(:,24).*g2(geometricdata(:,21),3))./geometricdata(:,18)];
%    
% g3d1 = [(geometricdata(:,8).*g3(geometricdata(:,5),1)+geometricdata(:,9).*g3(geometricdata(:,6),1)+geometricdata(:,10).*g3(geometricdata(:,7),1))./geometricdata(:,4)...
%         (geometricdata(:,8).*g3(geometricdata(:,5),2)+geometricdata(:,9).*g3(geometricdata(:,6),2)+geometricdata(:,10).*g3(geometricdata(:,7),2))./geometricdata(:,4)...
%         (geometricdata(:,8).*g3(geometricdata(:,5),3)+geometricdata(:,9).*g3(geometricdata(:,6),3)+geometricdata(:,10).*g3(geometricdata(:,7),3))./geometricdata(:,4)];
%     
% g3d2 = [(geometricdata(:,15).*g3(geometricdata(:,12),1)+geometricdata(:,16).*g3(geometricdata(:,13),1)+geometricdata(:,17).*g3(geometricdata(:,14),1))./geometricdata(:,11)...
%         (geometricdata(:,15).*g3(geometricdata(:,12),2)+geometricdata(:,16).*g3(geometricdata(:,13),2)+geometricdata(:,17).*g3(geometricdata(:,14),2))./geometricdata(:,11)...
%         (geometricdata(:,15).*g3(geometricdata(:,12),3)+geometricdata(:,16).*g3(geometricdata(:,13),3)+geometricdata(:,17).*g3(geometricdata(:,14),3))./geometricdata(:,11)];
% 
% g3d3 = [(geometricdata(:,22).*g3(geometricdata(:,19),1)+geometricdata(:,23).*g3(geometricdata(:,20),1)+geometricdata(:,24).*g3(geometricdata(:,21),1))./geometricdata(:,18)...
%         (geometricdata(:,22).*g3(geometricdata(:,19),2)+geometricdata(:,23).*g3(geometricdata(:,20),2)+geometricdata(:,24).*g3(geometricdata(:,21),2))./geometricdata(:,18)...
%         (geometricdata(:,22).*g3(geometricdata(:,19),3)+geometricdata(:,23).*g3(geometricdata(:,20),3)+geometricdata(:,24).*g3(geometricdata(:,21),3))./geometricdata(:,18)];
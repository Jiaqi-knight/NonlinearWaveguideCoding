function[tri1,tri2,tri3]=gettriangleadjacenttotriangle2D(tri,edges,triangles)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 24th, 2014
%    Last update: April 24th, 2014
%
%          Input: 
%         Output: 

%%
A = tri(1,7);
B = tri(1,8);
C = tri(1,9);

tri1 = 0;
tri2 = 0;
tri3 = 0;

[triBC1,triBC2] = gettriangleadjacenttoedge2D(edges(getedge_fromnodes2D(B,C,edges),:),triangles);

if triBC1~=tri(1,1)
    tri1 = triBC1;
elseif triBC2~=tri(1,1)
    tri1 = triBC2;
end

[triAB1,triAB2] = gettriangleadjacenttoedge2D(edges(getedge_fromnodes2D(A,B,edges),:),triangles);

if triAB1~=tri(1,1)
    tri2 = triAB1;
elseif triAB2~=tri(1,1)
    tri2 = triAB2;
end

[triAC1,triAC2] = gettriangleadjacenttoedge2D(edges(getedge_fromnodes2D(A,C,edges),:),triangles);

if triAC1~=tri(1,1)
    tri3 = triAC1;
elseif triAC2~=tri(1,1)
    tri3 = triAC2;
end

return
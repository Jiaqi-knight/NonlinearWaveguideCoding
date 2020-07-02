function[tri1,tri2]=gettriangleadjacenttoedge2D(edge,triangles)

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
A = edge(1,2);
B = edge(1,3);

tri1 = 0;
tri2 = 0;

i = 1;

while ~(tri1 && tri2) && i<=size(triangles,1)
    if ((triangles(i,7)==A && triangles(i,8)==B) || (triangles(i,7)==B && triangles(i,8)==A)) || ((triangles(i,7)==A && triangles(i,9)==B) || (triangles(i,7)==B && triangles(i,9)==A)) || ((triangles(i,8)==A && triangles(i,9)==B) || (triangles(i,8)==B && triangles(i,9)==A))
        if ~tri1
            tri1 = triangles(i,1);
        else
            tri2 = triangles(i,1);
        end
    end
    i = i + 1;
end

return
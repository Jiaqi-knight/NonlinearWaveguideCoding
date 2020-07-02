function[edge]=getedge_fromnodes2D(p1,p2,edges)

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
%          Input: 1 x 3 vector of point in plane (indices + coordinates)
%                 1 x 3 vector of point in plane (indices + coordinates)
%         Output: edge id

edge = 0;

i = 1;
search = 1;
while search && i<=size(edges,1)
    if (p1(1,1)==edges(i,2) && p2(1,1)==edges(i,3)) || (p1(1,1)==edges(i,3) && p2(1,1)==edges(i,2))
        edge = edges(i,1);
        search = 0;
    end
    i = i + 1;
end

return
function[l]=Rfunctioninterpolant(x,xval,index)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 14th, 2014
%    Last update: July 14th, 2014
%
%          Input: N x 1 vector x of interpolation nodes
%                 M x 1 vector xval of nodes for function evaluation
%                 index "index" of current node
%         Output: M x 1 vector l of interpolant evaluations
%%

N = size(x,1);
M = size(xval,1);

l = ones(M,1);

for i=1:N
    if i~=index
        l = l.*(xval-x(i,1))./(x(index,1)-x(i,1));
    end
end

return
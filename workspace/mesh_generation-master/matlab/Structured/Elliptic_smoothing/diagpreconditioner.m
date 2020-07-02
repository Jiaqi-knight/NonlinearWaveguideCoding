function[M,invM]=diagpreconditioner(A)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 19th, 2014
%    Last update: May 19th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

N = size(A,1);

M = zeros(N,N);
invM = zeros(N,N);

for i=1:N
    M(i,i) = A(i,i);
    invM(i,i) = 1/A(i,i);
end

return
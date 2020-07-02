function[M,invM]=SSORpreconditioner(A,omega)

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

D = zeros(N,N);
invD = zeros(N,N);
E = zeros(N,N);
F = zeros(N,N);

for i=1:N
    D(i,i) = A(i,i);
    invD(i,i) = 1/A(i,i);
    for j=1:N
        if j<i
            E(i,j) = - A(i,j);
            F(i,j) = 0;
        elseif j==i
            E(i,j) = 0;
            F(i,j) = 0;
        else
            E(i,j) = 0;
            F(i,j) = - A(i,j);
        end
    end
end

M = (D - omega*E)*invD*(D - omega*F);
invM = luinv(M);

return
function[M]=normaldistmthmoment(m)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 27th, 2014
%    Last update: May 27th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

n = floor(m/2);

M = zeros(n+1,3);

if 2*n==m
    for i = 0:n
        M(i+1,1) = binomialfactor(2*n,2*i)*normaldistcoeff(2*n-2*i);
        M(i+1,2) = 2*n-2*i;
        M(i+1,3) = 2*i;
    end
else
    for i = 0:n
        M(i+1,1) = binomialfactor(2*n+1,2*i+1)*normaldistcoeff((2*n+1)-(2*i+1));
        M(i+1,2) = (2*n+1)-(2*i+1);
        M(i+1,3) = 2*i+1;
    end
end

return
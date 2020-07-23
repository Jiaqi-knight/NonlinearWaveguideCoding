function[M]=evaluatenormaldistmthmoment(m,mu,sigma)

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

M = 0;

if 2*n==m
    for i = 0:n
        M = M + binomialfactor(2*n,2*i)*normaldistcoeff(2*n-2*i)*(sigma.^(2*n-2*i)).*(mu.^(2*i));
    end
else
    for i = 0:n
        M = M + binomialfactor(2*n+1,2*i+1)*normaldistcoeff((2*n+1)-(2*i+1))*(sigma.^((2*n+1)-(2*i+1))).*(mu.^(2*i+1));
    end
end

return
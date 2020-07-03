function[u]=computevelocity3D(N,Q,rho,f,scheme)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: June 24th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

u = zeros(N,3);

for i=1:Q
    u(:,1) = u(:,1) + f(:,i).*scheme(i,1);
    u(:,2) = u(:,2) + f(:,i).*scheme(i,2);
    u(:,3) = u(:,3) + f(:,i).*scheme(i,3);
end

u(:,1) = u(:,1)./rho;
u(:,2) = u(:,2)./rho;
u(:,3) = u(:,3)./rho;

return
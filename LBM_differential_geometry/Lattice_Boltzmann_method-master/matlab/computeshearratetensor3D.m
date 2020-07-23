function[S]=computeshearratetensor3D(N,Q,feq,f,rho,tau,scheme)

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

S = zeros(N,6);

for i=1:Q
    S(:,1) = S(:,1) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,1); % S11
    S(:,2) = S(:,2) + (f(:,i)-feq(:,i))*scheme(i,2)*scheme(i,2); % S22
    S(:,3) = S(:,3) + (f(:,i)-feq(:,i))*scheme(i,3)*scheme(i,3); % S33
    S(:,4) = S(:,4) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,2); % S12
    S(:,5) = S(:,5) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,3); % S13
    S(:,6) = S(:,6) + (f(:,i)-feq(:,i))*scheme(i,2)*scheme(i,3); % S23
end

S(:,1) = -(3./(2*tau(:,1))).*S(:,1)./rho(:,1);
S(:,2) = -(3./(2*tau(:,1))).*S(:,2)./rho(:,1);
S(:,3) = -(3./(2*tau(:,1))).*S(:,3)./rho(:,1);
S(:,4) = -(3./(2*tau(:,1))).*S(:,4)./rho(:,1);
S(:,5) = -(3./(2*tau(:,1))).*S(:,5)./rho(:,1);
S(:,6) = -(3./(2*tau(:,1))).*S(:,6)./rho(:,1);

return
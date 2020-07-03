function[sigma]=computestresstensor3D(N,Q,feq,f,tau,scheme)

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

sigma = zeros(N,6);

for i=1:Q
    sigma(:,1) = sigma(:,1) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,1); % sigma11
    sigma(:,2) = sigma(:,2) + (f(:,i)-feq(:,i))*scheme(i,2)*scheme(i,2); % sigma22
    sigma(:,3) = sigma(:,3) + (f(:,i)-feq(:,i))*scheme(i,3)*scheme(i,3); % sigma33
    sigma(:,4) = sigma(:,4) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,2); % sigma12
    sigma(:,5) = sigma(:,5) + (f(:,i)-feq(:,i))*scheme(i,1)*scheme(i,3); % sigma13
    sigma(:,6) = sigma(:,6) + (f(:,i)-feq(:,i))*scheme(i,2)*scheme(i,3); % sigma23
end

sigma(:,1) = -(1-0.5./tau(:,1)).*sigma(:,1);
sigma(:,2) = -(1-0.5./tau(:,1)).*sigma(:,2);
sigma(:,3) = -(1-0.5./tau(:,1)).*sigma(:,3);
sigma(:,4) = -(1-0.5./tau(:,1)).*sigma(:,4);
sigma(:,5) = -(1-0.5./tau(:,1)).*sigma(:,5);
sigma(:,6) = -(1-0.5./tau(:,1)).*sigma(:,6);

return
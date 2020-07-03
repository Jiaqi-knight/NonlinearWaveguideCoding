function[f]=collide3D(Q,fold,feq,Flambda,beta,dt)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: June 3rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

for q=1:Q
    f(:,q) = (1-beta(:,1)).*fold(:,q) + beta(:,1).*feq(:,q) + dt*Flambda(:,q);
end

return
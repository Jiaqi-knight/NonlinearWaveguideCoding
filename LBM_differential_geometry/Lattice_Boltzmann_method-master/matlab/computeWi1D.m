function[W]=computeWi1D(Q,velocities,T0)

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

if size(velocities,1)==1
    velocities = velocities';
end

u = 0;

T0sqrt = sqrt(T0);

M = zeros(Q,1);
vQ = zeros(Q,Q);

for i=0:Q-1
    M(i+1) = evaluatenormaldistmthmoment(i,u,T0sqrt);
    vQ(i+1,:) = velocities'.^i;
end

W = vQ\M;

return
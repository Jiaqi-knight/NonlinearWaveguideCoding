function[T]=T0closurerelation1D(Q,allvelocities)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 26th, 2014
%    Last update: May 26th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

velocities = allvelocities(2:end);

b = ones(0.5*(Q-1)+1,1);

for i=1:0.5*(Q-1)
    b(i+1) = (2*i+1)*b(i);
end

M = zeros(size(velocities,1),0.5*(Q-1)+1);

for i=0:0.5*(Q-1)
    M(:,i+1) = velocities.^(2*i+1);
end

a = M(:,1:end-1)\M(:,end);

a = [a;1];

c = a.*b;

c(1:end-1) = -c(1:end-1);

coeff = zeros(size(b,1),1);

for i=1:size(coeff)
    coeff(i,1) = c(size(c,1)-(i-1),1);
end

Troots = roots(coeff);

T = [];

for i=1:size(Troots,1)
    if isreal(Troots(i,1))
        T = [T;Troots(i,1)];
    end
end

return
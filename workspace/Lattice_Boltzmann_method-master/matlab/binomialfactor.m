function[fac]=binomialfactor(n,k)

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

fac = factorial(n)/(factorial(k)*factorial(n-k));

return
function[rho]=computerho3D(f)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: June 27th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

rho = sum(f,2);

return
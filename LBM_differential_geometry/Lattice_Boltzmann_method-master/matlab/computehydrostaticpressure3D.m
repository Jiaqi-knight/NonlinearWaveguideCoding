function[p]=computehydrostaticpressure3D(rho,cssq)

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

p = cssq.*rho;

return
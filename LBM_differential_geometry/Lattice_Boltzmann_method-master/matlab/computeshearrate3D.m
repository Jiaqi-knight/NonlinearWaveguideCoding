function[gamma]=computeshearrate3D(S)

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

S11 = S(:,1);
S22 = S(:,2);
S33 = S(:,3);
S12 = S(:,4);
S13 = S(:,5);
S23 = S(:,6);

gamma = sqrt(2.*(S11.*S11+S22.*S22+S33.*S33+2.*S12.*S12+2.*S13.*S13+2.*S23.*S23));

end
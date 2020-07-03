function[alpha]=normaldistcoeff(exp)

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

switch exp
    case 0
        alpha = 1;
    case 2
        alpha = 1;
    case 4
        alpha = 3;
    case 6
        alpha = 15;
    case 8
        alpha = 105;
    case 10
        alpha = 945;
    case 12
        alpha = 10395;
    case 14
        alpha = 135135;
    case 16
        alpha = 2027025;
    case 18
        alpha = 34459425;
    case 20
        alpha = 654729075;
    case 22
        alpha = 13749310575;
    case 24
        alpha = 316234143225;
    case 26
        alpha = 7905853580625;
end

return
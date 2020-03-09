clear all

%-------
% FSELIB
%
% diffusion matrix for a quadratic element
% with arbitrary interior node
%-------

beta = 0.8;

A(1,1) = 1/(1+beta)^2 * (4+3*(1+beta)^2 );
A(1,2) = -8/( (1+beta)*(1-beta^2) );
A(2,1) = A(1,2);
A(1,3) = 1/(1-beta^2) * (4-3*(1-beta^2));
A(3,1) = A(1,3);
A(2,2) = 16/(1-beta^2)^2;
A(2,3) = -8/((1-beta)*(1-beta^2));
A(3,2) = A(2,3);
A(3,3) = 1/(1-beta)^2 * (4 + 3*(1-beta)^2 );

A
det(A)

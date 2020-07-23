function [al, be, ga ] = elm6_abc ...
...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6)

%===============================
% Compute the (xi, eta) mapping coefficients
% alpha, beta, gamma
% for a six-node triangle
% in the plane
%===============================

D42 = sqrt( (x4-x2)^2 + (y4-y2)^2 );
D41 = sqrt( (x4-x1)^2 + (y4-y1)^2 );
D63 = sqrt( (x6-x3)^2 + (y6-y3)^2 );
D61 = sqrt( (x6-x1)^2 + (y6-y1)^2 );
D52 = sqrt( (x5-x2)^2 + (y5-y2)^2 );
D53 = sqrt( (x5-x3)^2 + (y5-y3)^2 );

al = 1.0/(1.0+D42/D41);
be = 1.0/(1.0+D63/D61);
ga = 1.0/(1.0+D52/D53);

%-----
% done
%-----

return;

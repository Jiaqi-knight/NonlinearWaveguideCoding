function [al, be, ga] = elm6_2d_abc ...
...
    (x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4 ...
    ,x5,y5,z5, x6,y6,z6)

%===============================
% Compute the (xi, eta) mapping
% coefficients alpha, beta, gamma
% for a six-node triangle
% in space
%===============================

D42 = sqrt( (x4-x2)^2 + (y4-y2)^2 + (z4-z2)^2 );
D41 = sqrt( (x4-x1)^2 + (y4-y1)^2 + (z4-z1)^2 );
D63 = sqrt( (x6-x3)^2 + (y6-y3)^2 + (z6-z3)^2 );
D61 = sqrt( (x6-x1)^2 + (y6-y1)^2 + (z6-z1)^2 );
D52 = sqrt( (x5-x2)^2 + (y5-y2)^2 + (z5-z2)^2 );
D53 = sqrt( (x5-x3)^2 + (y5-y3)^2 + (z5-z3)^2 );

al = 1.0/(1.0+D42/D41);
be = 1.0/(1.0+D63/D61);
ga = 1.0/(1.0+D52/D53);

%-----
% done
%-----

return;

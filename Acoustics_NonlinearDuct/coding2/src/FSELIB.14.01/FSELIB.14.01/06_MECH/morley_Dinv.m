function [Dinv] = morley_Dinv (x1,y1, x2,y2, x3,y3)

%====================================
% Compute the inverse of the matrix D
% for the Morley element
%====================================

%----------------------------------
% mid-side nodes and normal vectors
%----------------------------------

x4=0.5*(x1+x2); y4=0.5*(y1+y2);
x5=0.5*(x2+x3); y5=0.5*(y2+y3);
x6=0.5*(x3+x1); y6=0.5*(y3+y1);

d12 = sqrt((x2-x1)^2+(y2-y1)^2); nx4=-(y2-y1)/d12; ny4=(x2-x1)/d12;
d23 = sqrt((x3-x2)^2+(y3-y2)^2); nx5=-(y3-y2)/d23; ny5=(x3-x2)/d23;
d31 = sqrt((x1-x3)^2+(y1-y3)^2); nx6=-(y1-y3)/d31; ny6=(x1-x3)/d31;

%--------
% define the coefficient matrix
%--------

D = [ ...
1  x1  y1  x1^2  x1*y1  y1^2 ; ...
1  x2  y2  x2^2  x2*y2  y2^2 ; ...
1  x3  y3  x3^2  x3*y3  y3^2 ; ...
0 nx4 ny4 2*x4*nx4 y4*nx4+x4*ny4 2*y4*ny4 ; ...
0 nx5 ny5 2*x5*nx5 y5*nx5+x5*ny5 2*y5*ny5 ; ...
0 nx6 ny6 2*x6*nx6 y6*nx6+x6*ny6 2*y6*ny6 
];

%--------------------
% compute the inverse
%--------------------

Dinv = inv(D);

%-----
% done
%-----

return

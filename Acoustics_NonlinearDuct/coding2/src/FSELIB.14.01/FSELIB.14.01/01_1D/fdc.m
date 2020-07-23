%===============================
% FSELIB
%
% finite difference coefficients
%===============================

%============
% differentiation vector for (1.5.45)
% computed by the method of undetermined coefficients
%============

A= [ 1    1   1      1     1;
    -1    0   1      2     3;
    1/2   0  1/2    4/2   9/2;
   -1/6   0  1/6    8/6  27/6; 
    1/24  0  1/24  16/24 81/24];

b = [0 0 1 0 0];

sol1 = b/A';
sol1 = 4.0*sol1

[11 -20 6 4 -1]/3

%============
% differentiation vector for the third of (1.3.48)
% computed by the method of undetermined coefficients
%============

B= [ 1      1      1     1     1;
    -3     -2     -1     0     1;
    9/2    4/2    1/2    0   1/2;
  -27/6   -8/6   -1/6    0   1/6;
   81/24  16/24   1/24   0   1/24];

b=[0 0 1 0 0];

sol2 = b/B';
sol2 = 4.0*sol2

[-1 4 6 -20 11]/3

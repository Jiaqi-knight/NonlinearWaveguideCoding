clear all
close all

%============================================
% Driver for the
% 45-p Gaussian quadrature over a tetrahedron
%============================================

[xiq,etq,ztq,weq] = gauss_tetra45

sum = 0.0;
for i=1:45
 sum = sum+weq(i);
end

sum


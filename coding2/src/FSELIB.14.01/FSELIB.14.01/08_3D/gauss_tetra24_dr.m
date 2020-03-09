clear all
close all

%============================================
% Driver for the
% 24-p Gaussian quadrature over a tetrahedron
%============================================

[xiq,etq,ztq,weq] = gauss_tetra24;

sum = 0.0;
for i=1:24
 sum =sum+weq(i);
end

sum


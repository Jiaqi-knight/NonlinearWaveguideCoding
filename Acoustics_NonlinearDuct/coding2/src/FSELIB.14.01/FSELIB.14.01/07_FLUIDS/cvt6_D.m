function [Dx, Dy] = cvt6_D ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ)

%========================================
% Computation of the matrices D^x and D^y
% required by code cvt6
%========================================

%---------------------------------
% compute the mapping coefficients
%---------------------------------

[al, be, ga] = elm6_abc ...
...
  (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6);

%-----------------------------
% read the triangle quadrature
%-----------------------------

[xi, eta, w] = gauss_trgl(NQ);

%----------------------------------
% initialize the matrices Dx and Dy
%----------------------------------

 for k=1:6
   Dx(k) = 0.0; Dy(k) = 0.0;
 end

%-----------------------
% perform the quadrature
%-----------------------

for i=1:NQ

[psi, gpsi, hs] = elm6_interp ...
...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
    ,al,be,ga, xi(i),eta(i));

 cf = 0.5*hs*w(i);

 for k=1:6
   Dx(k) = Dx(k) + gpsi(k,1)*cf;
   Dy(k) = Dy(k) + gpsi(k,2)*cf;
 end

end

%-----
% done
%-----

return;

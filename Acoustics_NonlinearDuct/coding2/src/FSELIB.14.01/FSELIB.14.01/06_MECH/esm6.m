function [esm_xx, esm_xy, esm_yy, arel] = esm6 ...
...
     (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ)

%============================================
% Computation of the element stiffness matrix
% for plane stress analysis
%===========================

%---------------------------------
% compute the mapping coefficients
%---------------------------------

[al, be, ga] = elm6_abc ...
...
  (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6);

%---------
% read the triangle quadrature
%---------

[xi, eta, w] = gauss_trgl(NQ);

%---------
% initialize the stiffness matrices
%---------

 for k=1:6
  for l=1:6
   esm_xx(k,l) = 0.0;
   esm_xy(k,l) = 0.0;
   esm_yy(k,l) = 0.0;
  end
 end

%---------
% perform the quadrature
%---------

arel = 0.0;  % element area

for i=1:NQ

 [psi, gpsi, hs] = elm6_interp ...
 ...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
    ,al,be,ga, xi(i),eta(i));

 cf = 0.5*hs*w(i);

 for k=1:6
  for l=1:6
   esm_xx(k,l) = esm_xx(k,l) + gpsi(k,1)*gpsi(l,1)*cf;
   esm_xy(k,l) = esm_xy(k,l) + gpsi(k,1)*gpsi(l,2)*cf;
   esm_yy(k,l) = esm_yy(k,l) + gpsi(k,2)*gpsi(l,2)*cf;
  end
 end

 arel = arel + cf;

end

% disp (arel)

%-----
% done
%-----

return;

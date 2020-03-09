function [edm, arel] = edm6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ)

%===============================================
% Evaluation of the element diffusion matrix for
% a 6-node triangle using the Gauss-triangle
% integration quadrature with NQ base points
%===============================================

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

%----------------------------------------
% initialize the element diffusion matrix
%----------------------------------------

 for k=1:6
  for l=1:6
   edm(k,l) = 0.0;
  end
 end

%-----------------------
% perform the quadrature
%-----------------------

arel = 0.0;  % element area (optional)

for i=1:NQ

 [psi, gpsi, hs] = elm6_interp ...
 ...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
    ,al,be,ga, xi(i),eta(i));

 cf = 0.5*hs*w(i);

 for k=1:6
  for l=1:6
   edm(k,l) = edm(k,l) + (gpsi(k,1)*gpsi(l,1)   ...
                       +  gpsi(k,2)*gpsi(l,2) )*cf;
  end
 end

 arel = arel + cf;

end

% disp (arel)

%-----
% done
%-----

return;

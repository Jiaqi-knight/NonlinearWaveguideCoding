function [edm, emm, arel] = edmm6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ)

%======================================================
% Evaluation of the element diffusion and mass matrices
% for a 6-node triangle
% using a Gauss-triangle integration quadrature
%======================================================

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

%---------------------------------
% initialize the element diffusion 
% and mass matrices
%---------------------------------

 for k=1:6
  for l=1:6
   edm(k,l) = 0.0;
   emm(k,l) = 0.0;
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
   emm(k,l) = emm(k,l) + psi(k)*psi(l)*cf;
  end
 end

 arel = arel + cf;

end

% disp (arel)

%-----
% done
%-----

return;

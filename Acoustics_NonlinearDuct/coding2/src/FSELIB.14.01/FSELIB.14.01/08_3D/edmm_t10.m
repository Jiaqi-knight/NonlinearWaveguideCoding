function [edm, emm, volel] = edmm_t10 ...
...
   (x1,y1,z1, x2,y2,z2, x3,y3,z3 ...
  , x4,y4,z4, x5,y5,z5, x6,y6,z6 ...
   ,x7,y7,z7, x8,y8,z8, x9,y9,z9, x10,y10,z10)

%======================================================
% Evaluation of the element diffusion and mass matrices
% for a 10-node tetrahedron,
% using a Gauss-tetrahedron integration quadrature
%======================================================

%-----------------------------
% read the triangle quadrature
%-----------------------------

NQ = 45;
[xi, eta, zeta, w] = gauss_tetra45;

NQ = 24;
[xi, eta, zeta, w] = gauss_tetra24;

%---------------------------------
% initialize the element diffusion 
% and mass matrices
%---------------------------------

 for k=1:10
  for l=1:10
   edm(k,l) = 0.0;
   emm(k,l) = 0.0;
  end
 end

volel = 0.0;  % element volume (optional)

%-----------------------
% perform the quadrature
%-----------------------

for i=1:NQ

[psi, gpsi, hv] = elmt10_interp ...
...
   (x1,y1,z1, x2,y2,z2, x3,y3,z3 ...
   ,x4,y4,z4, x5,y5,z5, x6,y6,z6 ...
   ,x7,y7,z7, x8,y8,z8, x9,y9,z9, x10,y10,z10 ...
   ,xi(i), eta(i), zeta(i) );

 cf = hv*w(i)/6.0;

 for k=1:10
  for l=1:10
   edm(k,l) = edm(k,l) + (gpsi(k,1)*gpsi(l,1) ...
                       +  gpsi(k,2)*gpsi(l,2) ... 
                       +  gpsi(k,3)*gpsi(l,3) )*cf;
   emm(k,l) = emm(k,l) + psi(k)*psi(l)*cf;
  end
 end

 volel = volel + cf;

end

% disp (volel)

%edm
%det(edm)
%pause

%-----
% done
%-----

return;

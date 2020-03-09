function [edm, eam, arel] = edam6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
   ,u1,v1, u2,v2, u3,v3, u4,v4, u5,v5, u6,v6 ...
   ,NQ)

%================================================
% Evaluates the element diffusion and advection
% matrices and element area for a 6-node triangle
%================================================

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
% initialize the element diffusion 
% and advection matrices
%---------

 for k=1:6
  for l=1:6
   edm(k,l) = 0.0;
   eam(k,l) = 0.0;
  end
 end

%-----------------------
% perform the quadrature
%-----------------------

arel = 0.0;  % element area

for i=1:NQ

[psi, gpsi, hs] = elm6_interp ...
...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
    ,al,be,ga, xi(i),eta(i));

% interpolate the velocity at the quadrature base points

 u = u1*psi(1)+u2*psi(2)+u3*psi(3) ...
   + u4*psi(4)+u5*psi(5)+u6*psi(6);

 v = v1*psi(1)+v2*psi(2)+v3*psi(3) ...
   + v4*psi(4)+v5*psi(5)+v6*psi(6);

 cf = 0.5*hs*w(i);

 for k=1:6
  for l=1:6

   edm(k,l) = edm(k,l) + (gpsi(k,1)*gpsi(l,1)   ...
                       +  gpsi(k,2)*gpsi(l,2) )*hs*w(i);

   prj = u*gpsi(l,1)+v*gpsi(l,2);

   eam(k,l) = eam(k,l) + psi(k) * prj *cf;

  end
 end

 arel = arel + cf;

end

% disp (arel)

%-----
% done
%-----

return;

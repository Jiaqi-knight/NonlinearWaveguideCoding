function [gradx, grady] = elm6_grad (x,y,f)

%===================================
% Evaluate the node element gradient
%===================================

%---------------------------------
% compute the mapping coefficients
%---------------------------------

[al, be, ga] = elm6_abc ...
...
  (x(1),y(1), x(2),y(2), x(3),y(3)...
  ,x(4),y(4), x(5),y(5), x(6),y(6) );

%----------------------------
% nodal (xi, eta) coordinates
%----------------------------

xi(1)=0.0; eta(1)=0.0;  xi(2)=1.0; eta(2)=0.0;
xi(3)=0.0; eta(3)=1.0;  xi(4)=al;  eta(4)=0.0;
xi(5)=be ; eta(5)=1-be; xi(6)=0.0; eta(6)=ga;

%---------------------------
% compute the nodal gradient
%---------------------------

for l=1:6  % run over the nodes

 gradx(l) = 0.0; grady(l) = 0.0; 

  [psi, gpsi, hs] = elm6_interp ...
  ...
   (x(1),y(1), x(2),y(2), x(3),y(3) ...
   ,x(4),y(4), x(5),y(5), x(6),y(6) ...
   ,al,be,ga, xi(l),eta(l));

 for k=1:6
   gradx(l) = gradx(l)+f(k)*gpsi(k,1);
   grady(l) = grady(l)+f(k)*gpsi(k,2);
 end

end

%-----
% done
%-----

return;

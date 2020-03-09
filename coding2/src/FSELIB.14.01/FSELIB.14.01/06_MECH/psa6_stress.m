function [sig_xx, sig_xy, sig_yy] = psa6_stress ...
...
     (x,y,u,v,E,nu);

%================================
% Evaluation of the node element
% stresses in plane stress analysis
%================================

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
% compute the nodal stresses
%---------------------------

nus = nu^2;

for l=1:6  % run over the nodes

 sig_xx(l) = 0.0; sig_xy(l) = 0.0; sig_yy(l) = 0.0;

  [psi, gpsi, hs] = elm6_interp ...
  ...
   (x(1),y(1), x(2),y(2), x(3),y(3), x(4),y(4), x(5),y(5), x(6),y(6) ...
   ,al,be,ga, xi(l),eta(l));

 for k=1:6
   sig_xx(l) = sig_xx(l)...
       +(   u(k)*gpsi(k,1)+nu*v(k)*gpsi(k,2) ) /(1-nus);
   sig_xy(l) = sig_xy(l)...
       +0.5*(u(k)*gpsi(k,2)+v(k)*gpsi(k,1) ) /(1+nu);
   sig_yy(l) = sig_yy(l)...
       +(nu*u(k)*gpsi(k,1)+   v(k)*gpsi(k,2) ) /(1-nus);
 end

 sig_xx(l) = E*sig_xx(l);
 sig_xy(l) = E*sig_xy(l);
 sig_yy(l) = E*sig_yy(l);
 
end

%-----
% done
%-----

return;

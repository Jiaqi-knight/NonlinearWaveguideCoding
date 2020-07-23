%===============================================
% script gauss_trgl_dr
%
% Integrate over the standard triangle in xi-eta
% using the Gauss-triangle integration quadrature
%===============================================

NQ = 6;

%-----------------------------
% read the triangle quadrature
%-----------------------------

[xi, eta, w] = gauss_trgl(NQ);

%-----------------------
% perform the quadrature
%-----------------------

integral = 0.0; 

for i=1:NQ
 f = xi(i); % integrand
 f = eta(i); % integrand
 f = xi(i)^2; % integrand
 f = eta(i)^2; % integrand
 f = xi(i)*eta(i); % integrand
 f = xi(i)^3; % integrand
 f = eta(i)^3; % integrand
 f = xi(i)^2*eta(i); % integrand
 f = xi(i)^4; % integrand
 f = xi(i)^3*eta(i); % integrand
 f = xi(i)*eta(i)^3; % integrand
 f = xi(i)^2*eta(i)^2; % integrand
 f = xi(i)^5; % integrand
 integral = integral + 0.5*f*w(i);
end

disp (integral)

%-----
% done
%-----

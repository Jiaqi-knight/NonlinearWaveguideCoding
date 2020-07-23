%A function to compute given a current density vector either the polynomial
%or entropic quasiequilbria.
function equilibria = LBMQuasiEquilibria(rho,u,methodChoice,nx)
equilibria = zeros(3,nx);
if ( methodChoice == 0) %Polynomial Equilibria
equilibria(1,:) = rho./6.*(1 - 3*u + 3*u.^2);
equilibria(2,:) = 2.*rho./3.*(1 - 3*u.^2./2);
equilibria(3,:) = rho./6.*(1 + 3*u + 3*u.^2);
else %Entropic Equilibria
equilibria(1,:) = rho./6.*( -3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
equilibria(2,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*u.^2));
equilibria(3,:) = rho./6.*( 3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
end

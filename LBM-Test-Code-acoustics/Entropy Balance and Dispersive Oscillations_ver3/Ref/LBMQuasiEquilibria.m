%A function to compute given a current density vector either the polynomial
%or entropic quasiequilbria.
function [equilibria,T] = LBMQuasiEquilibria(scheme,rho,u,methodChoice,nx,T)
equilibria = zeros(size(scheme,1),nx);
if ( methodChoice == 0) %Polynomial Equilibria
equilibria(1,:) = rho./6.*(1 + 3*u + 3*u.^2);    
equilibria(2,:) = rho./6.*(1 - 3*u + 3*u.^2);
equilibria(3,:) = 2.*rho./3.*(1 - 3*u.^2./2);
else %Entropic Equilibria
equilibria1(1,:) = rho./6.*( 3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
equilibria1(2,:) = rho./6.*( -3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
equilibria1(3,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*u.^2));

%让我们来测试一下另外一种数值的办法，并对其进行优化加速
%Frapolli-ELBM4CompressibleFlow
w=scheme(:,2).';
c=scheme(:,1).';
[equilibria,T] =entropyEquilibrium(nx,1,size(scheme,1),rho.',w,c,u.',T);
equilibria=equilibria.';
max(max(equilibria1-equilibria))


end

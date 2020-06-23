
function [x,y,z,ds,dt] =hyp1setup()
zmin = -5.0; zmax = 5.0;
nz = 50; dz = (zmax-zmin)/nz;
thetamin = 0.0; thetamax = 2.0*pi;
ntheta = 50; dtheta = (thetamax-thetamin)/ntheta;
[paramz,theta] = meshgrid(zmin:dz:zmax, ...
thetamin:dtheta:thetamax+dtheta);
% Apply the parameterization to obtain R^3 coordinates
r = sqrt(1 + paramz.^2);
x = r .* cos(theta);
y = r .* sin(theta);
z = paramz;
ds=dz;
dt=dtheta;


surf(x,y,z);
end
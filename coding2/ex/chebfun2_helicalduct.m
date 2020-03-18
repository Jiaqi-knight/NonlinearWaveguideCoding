clc;clear;close all
subfunction_path2='.\chebfun-master'
addpath(genpath(subfunction_path2));
u = chebfun2(@(u,v) u, [0 6*pi 0 2*pi]);
v = chebfun2(@(u,v) v, [0 6*pi 0 2*pi]);
x = 2*(1-exp(u/(6*pi))).*cos(u).*cos(v/2).^2;
y = 2*(-1+exp(u/(6*pi))).*sin(u).*cos(v/2).^2;
z = 1-exp(u/(3*pi))-sin(v)+exp(u/(6*pi)).*sin(v);
surf(x,y,z), camlight
view(160,10), axis equal, box on
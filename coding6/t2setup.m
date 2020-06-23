% ------------------------------------------------------------------------%
% ----------Toroidal Surface ---------------------------------------------%
%-------------------------------------------------------------------------%

function [x,y,z,ds,dt]=t2setup()
R=5;
r=1;
u=linspace(0,2*pi,100);
v=linspace(0,2*pi,100);
ds=2*pi/100;
dt=2*pi/100;

[u,v]=meshgrid(u,v);

x=(R+r.*cos(v)).*cos(u);
y=(R+r.*cos(v)).*sin(u);
z=r.*sin(v);

figure(1)

mesh(x,y,z);
view([-52,64])

h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
xlabel('X','fontSize',14);
ylabel('Y','fontSize',14);
zlabel('Z','fontsize',14);
title('Torus','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white'); 
axis equal
end
%-------------------------------------------------------------------------
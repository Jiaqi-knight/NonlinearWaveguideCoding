function [surf_x,surf_y,surf_z,ds,dt]=helicalsetup()
s = 0:1:40;
h=exp(linspace(0,1,length(s)));
kappa=(0.2)./h;tau=0.4./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);
sw=sqrt(kappa.^2+tau.^2).*s;
mastX= a.*sin(sw);
mastY = a.*cos(sw);
mastZ = b.*sw;
ds=1;
dt=2*pi/51
[t]=frenet(mastX(:),mastY(:),mastZ(:));
[surf_x,surf_y,surf_z,t,n]=tubeplot(mastX,mastY,mastZ,h,s,50);
figure
surf(surf_x,surf_y,surf_z)
end

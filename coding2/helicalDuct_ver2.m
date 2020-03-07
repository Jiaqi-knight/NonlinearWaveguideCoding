% James-3.58
%h(s)
clc
clear
close all
s = 0:0.01:4;
h=0.1*exp(linspace(0,1.5,length(s)));
kappa=(2/3)./h;tau=0.2./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);

sw=sqrt(kappa.^2+tau.^2).*s;


n=1;
for k=1:n
x(k,:) = a.*sin(sw+2*pi*k/n);
y(k,:) = a.*cos(sw+2*pi*k/n);
z(k,:) = b.*sw;
end
% figure(1)
%  
% xlabel('x'); ylabel('y'); title('Circula helix');
% axis equal

figure
for k=1:n
tubeplot(x(k,:),y(k,:),z(k,:),h,s,50);hold on;
end
plot3(x, y, z);
daspect([1,1,1]); camlight;
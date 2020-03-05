% James-3.58
clc
clear
close all
h=0.1;
kappa=1/35/h;tau=0.1/h;
a=kappa/(kappa^2+tau^2);
b=tau/(kappa^2+tau^2);
s = 0:0.01:5;
sw=sqrt(kappa^2+tau^2)*s;


n=8;
for k=1:n
x(k,:) = a*sin(sw+2*pi*k/n);
y(k,:) = a*cos(sw+2*pi*k/n);
z(k,:) = b*sw;
end
% figure(1)
% plot3(x, y, z); 
% xlabel('x'); ylabel('y'); title('Circula helix');
% axis equal

figure(2)
for k=1:n
tubeplot(x(k,:),y(k,:),z(k,:),1*h,s,50);hold on;
end
daspect([1,1,1]); camlight;
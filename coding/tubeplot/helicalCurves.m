% Section 7.1: Helical curves
%
% Last Update: 03.08.2003
% t=linspace(0,2*pi,20);
% tubeplot([cos(t);sin(t);0.2*(t)],0.1);
% daspect([1,1,1]); camlight;
% ---------------------------------------------------------------------------- %
% James-3.58
h=0.1;
kappa=2/3/h;tau=0.2/h;
a=kappa/(kappa^2+tau^2);
b=tau/(kappa^2+tau^2);
s = 0:0.001:1;
sw=sqrt(kappa^2+tau^2)*s;

x = a*sin(sw);
y = a*cos(sw);
z = b*sw;

figure(1)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Circula helix');
axis equal

tubeplot([x;y;z],h,51,0.01);
daspect([1,1,1]); camlight;
% % ---------------------------------------------------------------------------- %
% % 7.1.2: Elliptical helix
% a = 0.3;
% b = 1.0;
% c = 5.0;
% s = 0:0.01:10*pi;
% 
% x = a*sin(s);
% y = b*cos(s);
% z = s/(2*pi*c);
% 
% figure(2);
% plot3(x, y, z);
% xlabel('x'); ylabel('y'); title('Elliptical helix');
%     
% % ---------------------------------------------------------------------------- %
% % 7.1.3: Conical helix
% a = 0.5;
% c = 5.0;
% s = 0:0.01:10*pi;
% 
% x = (a*s/2*pi*c).*sin(s);
% y = (a*s/2*pi*c).*cos(s);
% z = s/(2*pi*c);
% 
% figure(3);
% plot3(x, y, z);
% xlabel('x'); ylabel('y'); title('Conical helix');
% 
% % ---------------------------------------------------------------------------- %
% % 7.1.4: Spherical helix
% c = 5.0;
% s = 0:0.01:10*pi;
% 
% x = sin(s/(2*c)).*cos(s);
% y = sin(s/(2*c)).*sin(s);
% z = cos(s/(2*c));
% 
% figure(4);
% plot3(x, y, z);
% xlabel('x'); ylabel('y'); title('Spherical helix');
% 
% 
% % ---------------------------------------------------------------------------- %
% % 7.1.5: n-Helix
% a = 0.3;
% c = 3.0;
% n = 3;
% s = 0:0.01:6*pi;
% 
% x1 = a*cos(s+2*pi*1/n);
% y1 = a*sin(s+2*pi*1/n);
% x2 = a*cos(s+2*pi*2/n);
% y2 = a*sin(s+2*pi*2/n);
% x3 = a*cos(s+2*pi*3/n);
% y3 = a*sin(s+2*pi*3/n);
% z = cos(s/(2*c));
% 
% figure(5);
% hold on;
% plot3(x1, y1, z);
% plot3(x2, y2, z, 'g-');
% plot3(x3, y3, z, 'r-');
% xlabel('x'); ylabel('y'); title('3-helix');
% hold off;
% view(3);

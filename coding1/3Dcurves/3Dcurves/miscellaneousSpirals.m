% Section 7.3: Miscellaneous spirals
%
% Last Update: 03.08.2003

% ---------------------------------------------------------------------------- %
% 7.3.1: SiCi spiral
% (Symbolic Math Toolbox required)
a = 0.5;
c = 20.0;
t = 0.2:0.01:20;

x = a*cosint(t);
y = a*sinint(t);
z = t/c;

figure(1)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('SiCi spiral');

% ---------------------------------------------------------------------------- %
% 7.3.3: Toroidal spiral
a = 0.2;
b = 0.8;
c = 20.0;
t = 0:0.01:2*pi;

x = (a*sin(c*t)+b).*cos(t);
y = (a*sin(c*t)+b).*sin(t);
z = a*cos(c*t);

figure(2)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Toroidal spiral');

% ---------------------------------------------------------------------------- %
% 7.3.4: Intersection of Sphere and Cylinder
t = 0:0.01:2*pi;

x = (1+cos(t))/2;
y = sin(t)/2
z = sin(t/2);

figure(3)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Intersection of Sphere and Cylinder');

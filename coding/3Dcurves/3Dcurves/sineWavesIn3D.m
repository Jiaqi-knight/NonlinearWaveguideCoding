% Section 7.2: Sine waves in three dimensions
%
% Last Update: 03.08.2003

% ---------------------------------------------------------------------------- %
% 7.2.1: Sine wave on Cylinder
a = 10.0;
b = 1.0;
c = 0.5;
t = 0:0.01:2*pi;

x = b*cos(t);
y = b*sin(t);
z = c*cos(a*t);

figure(1)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Sine wave on Cylinder');

% ---------------------------------------------------------------------------- %
% 7.2.2: Sine wave on Sphere
a = 10.0;
b = 1.0;
c = 0.3;
t = 0:0.01:2*pi;

x = (b^2 - c^2*cos(a*t).^2).^(1/2).*cos(t);
y = (b^2 - c^2*cos(a*t).^2).^(1/2).*sin(t);
z = c*cos(a*t);

figure(2)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Sine wave on Cylinder');

% ---------------------------------------------------------------------------- %
% 7.2.3: Sine wave on Hyperboloid of one sheet
a = 10.0;
b = 1.0;
c = 0.3;
t = 0:0.01:2*pi;

x = (b^2 + c^2*cos(a*t).^2).^(1/2).*cos(t);
y = (b^2 + c^2*cos(a*t).^2).^(1/2).*sin(t);
z = c*cos(a*t);

figure(3)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Sine wave on Hyperboloid of one sheet');

% ---------------------------------------------------------------------------- %
% 7.2.4: Sine wave on Cone
a = 10.0;
b = 0.5;
c = 0.4;
t = 0:0.01:2*pi;

x = b*(1+cos(a*t)).*cos(t);
y = b*(1+cos(a*t)).*sin(t);
z = c*(1+cos(a*t));

figure(4)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Sine wave on Cone');

% ---------------------------------------------------------------------------- %
% 7.2.5: Rotating sine wave
a = 3.0;
b = 1.0;
c = 1.0;
t = -2*pi:0.01:2*pi;

x = sin(a*t).*cos(b*t);
y = sin(a*t).*sin(b*t);
z = c*t/(2*pi);

figure(5)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Rotating sine wave 1');

a = 3.0;
b = 0.25;
c = 1.0;
t = -2*pi:0.01:2*pi;

x = sin(a*t).*cos(b*t);
y = sin(a*t).*sin(b*t);
z = c*t/(2*pi);

figure(6)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Rotating sine wave 2');

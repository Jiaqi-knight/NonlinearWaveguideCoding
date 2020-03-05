% Section 7.1: Helical curves
%
% Last Update: 03.08.2003

% ---------------------------------------------------------------------------- %
% 7.1.1: Circular helix (right helicoid)
a = 0.5;
c = 5.0;
t = 0:0.01:10*pi;

x = a*sin(t);
y = a*cos(t);
z = t/(2*pi*c);

figure(1)
plot3(x, y, z); 
xlabel('x'); ylabel('y'); title('Circula helix');

% ---------------------------------------------------------------------------- %
% 7.1.2: Elliptical helix
a = 0.3;
b = 1.0;
c = 5.0;
t = 0:0.01:10*pi;

x = a*sin(t);
y = b*cos(t);
z = t/(2*pi*c);

figure(2);
plot3(x, y, z);
xlabel('x'); ylabel('y'); title('Elliptical helix');
    
% ---------------------------------------------------------------------------- %
% 7.1.3: Conical helix
a = 0.5;
c = 5.0;
t = 0:0.01:10*pi;

x = (a*t/2*pi*c).*sin(t);
y = (a*t/2*pi*c).*cos(t);
z = t/(2*pi*c);

figure(3);
plot3(x, y, z);
xlabel('x'); ylabel('y'); title('Conical helix');

% ---------------------------------------------------------------------------- %
% 7.1.4: Spherical helix
c = 5.0;
t = 0:0.01:10*pi;

x = sin(t/(2*c)).*cos(t);
y = sin(t/(2*c)).*sin(t);
z = cos(t/(2*c));

figure(4);
plot3(x, y, z);
xlabel('x'); ylabel('y'); title('Spherical helix');


% ---------------------------------------------------------------------------- %
% 7.1.5: n-Helix
a = 0.3;
c = 3.0;
n = 3;
t = 0:0.01:6*pi;

x1 = a*cos(t+2*pi*1/n);
y1 = a*sin(t+2*pi*1/n);
x2 = a*cos(t+2*pi*2/n);
y2 = a*sin(t+2*pi*2/n);
x3 = a*cos(t+2*pi*3/n);
y3 = a*sin(t+2*pi*3/n);
z = cos(t/(2*c));

figure(5);
hold on;
plot3(x1, y1, z);
plot3(x2, y2, z, 'g-');
plot3(x3, y3, z, 'r-');
xlabel('x'); ylabel('y'); title('3-helix');
hold off;
view(3);

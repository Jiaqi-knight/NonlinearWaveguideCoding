clear all
close all

e11=1.0; e12=0.0; e13=0.0;
e21=0.0; e22=1.0; e23=0.0;
e31=0.0; e32=0.0; e33=1.0;

a11 = 0.5*(-e11+e21+e31);
a12 = 0.5*(-e12+e22+e32);
a13 = 0.5*(-e13+e23+e33);

a21 = 0.5*( e11-e21+e31);
a22 = 0.5*( e12-e22+e32);
a23 = 0.5*( e13-e23+e33);

a31 = 0.5*( e11+e21-e31);
a32 = 0.5*( e12+e22-e32);
a33 = 0.5*( e13+e23-e33);

figure(1)
colormap([0 0 0]);
axis equal
axis([-1 1 -1 1 -1 1]);

[x,y,z]=sphere(16);

radius = 0.25;
x = radius*x; y=radius*y; z=radius*z;

surf(x,y,z);
hold on;

u = ones(size(x));

%-------------

for i=-1:1
for j=-1:1
for k=-1:1

x1=x+(i*a11+j*a21+k*a31)*u;
y1=y+(i*a12+j*a22+k*a32)*u;
z1=z+(i*a13+j*a23+k*a33)*u;
axis([-1 1 -1 1 -1 1]);
 mesh(x1,y1,z1);
 surf(x1,y1,z1);

end
end
end


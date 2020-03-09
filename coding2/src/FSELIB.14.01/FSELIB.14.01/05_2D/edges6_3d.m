
clear all
close all

%=========================================
% visualize the edges of a 6-node triangle
%=========================================

figure(1)
hold on
axis equal
box on
set(gca,'fontsize',14)
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('z','fontsize',14);
view([77,10])

%============
for repeat=1:2
%============

%-------------
% one triangle
%-------------

if(repeat==1)

x1=-0.1; y1=-0.2; z1=-0.3;
x2= 0.93; y2= 0.10; z2=0.3;
x3= 0.05; y3= 1.15; z3=0.1;
x4= 0.8; y4=-0.1; z4=-0.1;
x5= 0.7;  y5= 0.5; z5= 0.3;
x6=-0.1; y6=0.1;  z6=0.6;
end

%-----------------
% another triangle
%-----------------

if(repeat==2)

x1=0.7;  y1= 1.50; z1=-0.8;
x2=0.05; y2= 1.15; z2=0.1;
x3=0.93; y3= 0.10; z3=0.3;
x4=0.5;  y4= 1.1;  z4=-0.1;
x5=0.7;  y5= 0.5;  z5= 0.3;
x6=0.83; y6= 0.5;  z6=-0.2;

end

%-----------------------
% relabel the nodes once
%-----------------------

% save
% x1s=x1; y1s=y1; x2s=x2; y2s=y2; x3s=x3; y3s=y3;
% x4s=x4; y4s=y4; x5s=x5; y5s=y5; x6s=x6; y6s=y6;

% relabel
% x1=x2s; y1=y2s; x2=x3s; y2=y3s; x3=x1s; y3=y1s;
% x4=x5s; y4=y5s; x5=x6s; y5=y6s; x6=x4s; y6=y4s;

%------------------------
% relabel the nodes twice
%------------------------

% save
% x1s=x1; y1s=y1; x2s=x2; y2s=y2; x3s=x3; y3s=y3;
% x4s=x4; y4s=y4; x5s=x5; y5s=y5; x6s=x6; y6s=y6;

% relabel
% x1=x2s; y1=y2s; x2=x3s; y2=y3s; x3=x1s; y3=y1s;
% x4=x5s; y4=y5s; x5=x6s; y5=y6s; x6=x4s; y6=y4s;

%---------------------------------------------
% compute the coefficients: alpha, beta, gamma
%---------------------------------------------

[al,be,ga] = elm6_3d_abc ...
...
    (x1,y1,z1, x2,y2,z2, x3,y3,z3 ...
    ,x4,y4,z4, x5,y5,z5, x6,y6,z6);

%al = 0.5;
%be = 0.5;
%ga = 0.5;

%--------
% prepare
%--------

np = 32;
step = 1.0/np;

% al=0.5; be=0.5; ga=0.5;

alc = 1.0-al;
bec = 1.0-be;
gac = 1.0-ga;

%---------------
% plot the nodes
%---------------

plot3 (x1,y1,z1,'ko');
plot3 (x2,y2,z2,'ko');
plot3 (x3,y3,z3,'ko');
plot3 (x4,y4,z4,'ko');
plot3 (x5,y5,z5,'ko');
plot3 (x6,y6,z6,'ko');

% patch ([x1,x4,x2,x5,x3,x6,x1], ...
%       [y1,y4,y2,y5,y3,y6,y1], ...
%       [z1,z4,z2,z5,z3,z6,z1], ...
%       [z1,z4,z2,z5,z3,z6,z1])

%-----------
% edge 1-4-2
%-----------

for i=1:np+1
 xi=(i-1.0)*step;
 psi2 = xi*(xi-al)/alc ;
 psi4 = xi*(1.0-xi)/(al*alc) ;
 psi1 = 1.0-psi2-psi4;
 x(i) = x1*psi1 + x2*psi2 + x4*psi4;
 y(i) = y1*psi1 + y2*psi2 + y4*psi4;
 z(i) = z1*psi1 + z2*psi2 + z4*psi4;
end

%plot3 (x,y,'--g');
plot3 (x,y,z,'k');

%-----------
% edge 1-6-3
%-----------

for i=1:np+1
 eta=(i-1.0)*step;
 psi3 = eta*(eta-be)/bec ;
 psi6 = eta*(1.0-eta)/(be*bec) ;
 psi1 = 1.0-psi3-psi6;
 x(i) = x1*psi1 + x3*psi3 + x6*psi6;
 y(i) = y1*psi1 + y3*psi3 + y6*psi6;
 z(i) = z1*psi1 + z3*psi3 + z6*psi6;
end

%plot3 (x,y,z,'--r');
plot3 (x,y,z,'k');

%-----------
% edge 2-5-3
%-----------

for i=1:np+1
 xi=(i-1.0)*step;
 eta = 1.0-xi;
 psi2 = xi*(xi-al+eta*(al-ga)/gac)/alc;
 psi5 = xi*eta/(ga*gac);
 psi3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec;
 x(i) = x2*psi2 + x3*psi3 + x5*psi5;
 y(i) = y2*psi2 + y3*psi3 + y5*psi5;
 z(i) = z2*psi2 + z3*psi3 + z5*psi5;
end

%plot3 (x,y,'--c');
plot3 (x,y,z,'k');

%======
end  % of repeat
%======

%-----
% done
%-----

%=======
% morley
%
% prepare graphs of the element interpolation functions
% for the Morley element
%=======

clear all
close all

%--------------------------------
% inquire the mode to be graphed
%-------------------------------

mode = input('Please enter the mode to plot: ')

%-------
% define the vertices arbitrarily
%-------

x1= 0.0; y1= 0.0; x2= 1.0; y2= 0.0; x3= 0.0; y3= 1.0;
%x1= 1.0; y1= 0.0; x2= 2.0; y2= 1.5; x3= 0.0; y3= 1.0;

%-------------------------------------------
% compute the inverse matrix of the D matrix
%-------------------------------------------

Dinv = morley_Dinv (x1,y1, x2,y2, x3,y3);

%-----------
% prepare
%-----------

figure(1)
hold on
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
view([58, 26])
%axis equal
box on

%---------------------------------------------
% define the plotting grid in the xi-eta plane
%---------------------------------------------

N=32; M=32;  % arbitrary

Dxi = 1.0/N; Deta = 1.0/M;

for i=1:N+1 
  xi(i)=Dxi*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end

for j=1:M+1
  eta(j)=Deta*(j-1.0);
% fprintf (1,'%6.2f',Y(j))
end

%----------------------
% horizontal grid lines
%----------------------

for i=1:N
for j=1:M+2-i
 x = x1 + (x2-x1)*xi(i) + (x3-x1)*eta(j);
 y = y1 + (y2-y1)*xi(i) + (y3-y1)*eta(j);
 monvec = [1  x y  x^2 x*y y^2];
 psi = monvec*Dinv;
 xp(j) = x; yp(j)=y; zp(j) = psi(mode);
end
plot3(xp,yp,zp,'k');
end

%----------------------
% vertical grid lines
%----------------------

for j=1:M
for i=1:N+2-j
 x = x1 + (x2-x1)*xi(i) + (x3-x1)*eta(j);
 y = y1 + (y2-y1)*xi(i) + (y3-y1)*eta(j);
 monvec = [1  x y  x^2 x*y y^2];
 psi = monvec*Dinv;
 xp(i) = x; yp(i)=y; zp(i) = psi(mode);
end
plot3(xp,yp,zp,'k');
end

%------------
% plot the perimeter of the triangle
%-----------

xx(1)=x1;yy(1)=y1;zz(1)=0.0;
xx(2)=x2;yy(2)=y2;zz(2)=0.0;
xx(3)=x3;yy(3)=y3;zz(3)=0.0;
xx(4)=x1;yy(4)=y1;zz(4)=0.0;

plot3(xx,yy,zz,'k');

%-----
% done
%-----

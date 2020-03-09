%=======
% HCT
%
% prepare graphs of the interpolation functions
% for the HCT element
%=======

clear all
close all

%--------------------------------
% inquire the mode to be graphed
%-------------------------------

mode = input('Please enter the mode to plot (1-12): ')

dof = zeros(1,12);
dof(mode)=1.0;

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

%-------
% define the vertices arbitrarily
%-------

x1= 0.0; y1= 0.0; x2= 1.0; y2= 0.0; x3= 0.5; y3= sqrt(3.0)/2.0;
x1= -0.5; y1= 0.0; x2= 0.5; y2= 0.0; x3= 0.0; y3= sqrt(3.0)/2.0;
%x1= 1.0; y1= 0.0; x2= 2.0; y2= 1.5; x3= 0.0; y3= 1.0;
%x1= 1.0; y1= 1.0; x2= 2.0; y2= 1.0; x3= 1.0; y3= 2.0;
x1= 0.0; y1= 0.0; x2= 1.0; y2= 0.0; x3= 0.0; y3= 1.0;

%-----------------------------------------
% Generate the modal polynomial coefficients
% for the HCTF element
%-----------------------------------------

[a] = HCT_sys (x1,y1, x2,y2, x3,y3, dof);

%----------------
% proceed to plot
%----------------

x7=(x1+x2+x3)/3.0;
y7=(y1+y2+y3)/3.0;

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

%======
for pass=1:3
%======

%----------------------
% horizontal grid lines
%----------------------

for i=1:N
for j=1:M+2-i
 if(pass==1)
   x = x1 + (x2-x1)*xi(i) + (x7-x1)*eta(j);
   y = y1 + (y2-y1)*xi(i) + (y7-y1)*eta(j);
   zp(j) = a(1)+a(2)*x+a(3)*y+a(4)*x^2+a(5)*x*y+a(6)*y^2 ...
             +a(7)*x^3+a(8)*x^2*y+a(9)*x*y^2+a(10)*y^3;
 elseif(pass==2)
   x = x2 + (x3-x2)*xi(i) + (x7-x2)*eta(j);
   y = y2 + (y3-y2)*xi(i) + (y7-y2)*eta(j);
   zp(j) = a(11)+a(12)*x+a(13)*y+a(14)*x^2+a(15)*x*y+a(16)*y^2 ...
                +a(17)*x^3+a(18)*x^2*y+a(19)*x*y^2+a(20)*y^3;
 else
   x = x3 + (x1-x3)*xi(i) + (x7-x3)*eta(j);
   y = y3 + (y1-y3)*xi(i) + (y7-y3)*eta(j);
   zp(j) = a(21)+a(22)*x+a(23)*y+a(24)*x^2+a(25)*x*y+a(26)*y^2 ...
                +a(27)*x^3+a(28)*x^2*y+a(29)*x*y^2+a(30)*y^3;
 end
 xp(j)=x; yp(j)=y;
end
plot3(xp,yp,zp,'k'); 
end

%----------------------
% vertical grid lines
%----------------------

for j=1:M
for i=1:N+2-j
 if(pass==1)
   x = x1 + (x2-x1)*xi(i) + (x7-x1)*eta(j);
   y = y1 + (y2-y1)*xi(i) + (y7-y1)*eta(j);
   zp(i) = a(1)+a(2)*x+a(3)*y+a(4)*x^2+a(5)*x*y+a(6)*y^2 ...
               +a(7)*x^3+a(8)*x^2*y+a(9)*x*y^2+a(10)*y^3;
 elseif(pass==2)
   x = x2 + (x3-x2)*xi(i) + (x7-x2)*eta(j);
   y = y2 + (y3-y2)*xi(i) + (y7-y2)*eta(j);
   zp(i) = a(11)+a(12)*x+a(13)*y+a(14)*x^2+a(15)*x*y+a(16)*y^2 ...
                +a(17)*x^3+a(18)*x^2*y+a(19)*x*y^2+a(20)*y^3;
 else
   x = x3 + (x1-x3)*xi(i) + (x7-x3)*eta(j);
   y = y3 + (y1-y3)*xi(i) + (y7-y3)*eta(j);
   zp(i) = a(21)+a(22)*x+a(23)*y+a(24)*x^2+a(25)*x*y+a(26)*y^2 ...
                +a(27)*x^3+a(28)*x^2*y+a(29)*x*y^2+a(30)*y^3;
 end
 xp(i)=x; yp(i)=y;
end
plot3(xp,yp,zp,'k');
end

%======
end
%======

%------------
% plot the perimeter of the triangle
%-----------

xx(1)=x1;yy(1)=y1;zz(1)=0.0;
xx(2)=x2;yy(2)=y2;zz(2)=0.0;
xx(3)=x3;yy(3)=y3;zz(3)=0.0;
xx(4)=x1;yy(4)=y1;zz(4)=0.0;

plot3(xx,yy,zz,'k');

%-----------------------
% plot the interior node
%-----------------------

plot3(x7,y7,0.0,'ko');

%-------------------------
% bending stiffness matrix (optional)
%-------------------------

ido =0;

if(ido==1)

NQ=4;
[ebsm, rhs, arel] = HCT_ebsm (x1,y1, x2,y2, x3,y3, NQ);

end

%-----
% done
%-----


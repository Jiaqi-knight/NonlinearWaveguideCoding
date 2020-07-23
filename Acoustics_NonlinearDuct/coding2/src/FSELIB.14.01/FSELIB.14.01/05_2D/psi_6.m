
clear all
close all

%================================
% psi_6
%
% prepare graphs of the element
% interpolation functions psi
% for a quadratic, 6-node triangle
%=================================

node = 1;
node = 4;

%---
% prepare
%---

figure(1)
axis equal
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);

if(node==1)
 zlabel('\psi_1','fontsize',14);
 axis([0 1 0 1 -0.2 1])
elseif(node==4)
 zlabel('\psi_4','fontsize',14);
 axis([0 1 0 1 0 1])
end

set(gca,'fontsize',14)
box on
view([60,30])

N=32; M=32;   % plotting grid

Dx=1.0/N;
Dy=1.0/M;

for i=1:N+1 
  x(i)=Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end

for j=1:M+1
  y(j)=Dy*(j-1.0);
% fprintf (1,'%6.2f',Y(j))
end

%---
% vertical lines
%---

for i=1:N
 for j=1:M+2-i
  xp(j) = x(i);
  yp(j) = y(j);
  zeta = 1.0-xp(j)-yp(j);
  if(node==1)
   zp(j)= zeta*(2.0*zeta-1.0);   % 6-node vertex
  elseif(node==4)
   zp(j)= 4.0*xp(j)*zeta;        % 6-node mid
  end
 end
plot3(xp,yp,zp,'k');
end

clear xp yp zp

%---
% horizontal lines
%---

for j=1:M
 for i=1:N+2-j
  yp(i) = y(j);
  xp(i) = x(i);
  zeta = 1.0-xp(i)-yp(i);
  if(node==1)
   zp(i)= zeta*(2.0*zeta-1.0);  % 6-node vertex
  elseif(node==4)
   zp(i)= 4.0*xp(i)*zeta;       % 6-node mid
  end
 end
plot3(xp,yp,zp,'k')
end

%---
% plot the triangle
%---

xx(1)=0.0; yy(1)=0.0; zz(1)=0.0;
xx(2)=1.0; yy(2)=0.0; zz(2)=0.0;
xx(3)=0.0; yy(3)=1.0; zz(3)=0.0;
xx(4)=0.0; yy(4)=0.0; zz(4)=0.0;

plot3(xx,yy,zz,'k')
                                                                                
%-----
% done
%-----

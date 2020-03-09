clear all
close all

%=========================
% generate a brick element
% defined by 8 nodes
%=========================

%===
% prepare
%===

figure(1)
hold on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)
set(gca,'fontsize',14)
axis('equal')
box on

%-----
% vertex coordinates (typical)
%-----

x(1) = 0.0; y(1)= 0.0; z(1)= 0.00;
x(2) = 0.9; y(2)= 0.2; z(2)= 0.02;
x(3) = 0.1; y(3)= 1.1; z(3)=-0.05;
x(4) = 0.8; y(4)= 1.0; z(4)= 0.02;
x(5) = 0.1; y(5)= 0.1; z(5)= 1.01;
x(6) = 0.8; y(6)= 0.6; z(6)= 1.10;
x(7) = 0.2; y(7)= 0.9; z(7)= 0.90;
x(8) = 0.7; y(8)= 0.8; z(8)= 0.80;

%-----
% face node labels
%-----

s(1,1)=1; s(1,2)=2; s(1,3)=4; s(1,4)=3;
s(2,1)=1; s(2,2)=3; s(2,3)=7; s(2,4)=5;
s(3,1)=6; s(3,2)=5; s(3,3)=7; s(3,4)=8;
s(4,1)=2; s(4,2)=6; s(4,3)=8; s(4,4)=4;
s(5,1)=3; s(5,2)=4; s(5,3)=8; s(5,4)=7;
s(6,1)=2; s(6,2)=6; s(6,3)=5; s(6,4)=1;

%-----
% paint the faces
%-----

for i=1:6
 for j=1:4
  xp(j)=x(s(i,j)); yp(j)=y(s(i,j)); zp(j)=z(s(i,j));
 end
   if(i==1) patch(xp,yp,zp,'r'); end
   if(i==2) patch(xp,yp,zp,'b'); end
   if(i==3) patch(xp,yp,zp,'y'); end
   if(i==4) patch(xp,yp,zp,'g'); end
   if(i==5) patch(xp,yp,zp,'r'); end
   if(i==6) patch(xp,yp,zp,'b'); end
%  plot3(xp,yp,zp); hold on;
end

%-----
% done
%-----

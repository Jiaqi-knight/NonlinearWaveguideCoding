clear all
close all

%================================================
% Code nodes_tetra_lob
%
% Lobatto node distribution in the tetrahedron
% corresponding to a complete mth order expansion
%
% The nodes over the faces are placed at the
% Lobatto triangle set
%================================================

%----
% prepare
%----

figure(1)
hold on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)
set(gca,'fontsize',14)
axis('equal')
axis([0 1 0 1 0 1])
box on
view(56,22)

%-------------
% parameters
%-------------

m = 3
m = 4

%-------------
% lobatto grid
%-------------

[t, W] = lobatto(m-1);    % W will not be needed
                                                                                
v(1) = 0.0;
for i=2:m
 v(i) = 0.5*(1.0+t(i-1));
end
v(m+1) = 1.0;

%--------------------
% uniform master grid
%--------------------

%for i=1:m+1
%  v(i) = (i-1.0)/m;
%end

%---------------------------
% deploy and count the nodes
% on the faces
% uncomment for a regular tetrahedron
%---------------------------

% face z=0:

k=1;
for i=1:m+1;
 for j=1:m+2-i
     l = m+4-i-j-k;
     xp = (1.0+2.0*v(i)-    v(j)-    v(k)-v(l))/3.0;
     yp = (1.0-    v(i)+2.0*v(j)-    v(k)-v(l))/3.0;
     zp = 0;
%        xp=xp+0.5*(yp+zp); yp=0.5*sqrt(3)*(yp+zp/3.0);
%        zp = sqrt(2/3)*zp;
     plot3(xp,yp,zp,'ko');
  end
end

% face x=0:

k=1;
for i=1:m+1;
 for j=1:m+2-i
     l = m+4-i-j-k;
     xp = 0;
     yp = (1.0+2.0*v(i)-    v(j)-    v(k)-v(l))/3.0;
     zp = (1.0-    v(i)+2.0*v(j)-    v(k)-v(l))/3.0;
%        xp=xp+0.5*(yp+zp); yp=0.5*sqrt(3)*(yp+zp/3.0);
%        zp = sqrt(2/3)*zp;
     plot3(xp,yp,zp,'x'); hold on;
  end
end

% face y=0:

k=1;
for i=1:m+1;
 for j=1:m+2-i
     l = m+4-i-j-k;
     xp = (1.0+2.0*v(i)-    v(j)-    v(k)-v(l))/3.0;
     yp = 0;
     zp = (1.0-    v(i)+2.0*v(j)-    v(k)-v(l))/3.0;
%        xp=xp+0.5*(yp+zp); yp=0.5*sqrt(3)*(yp+zp/3.0);
%        zp = sqrt(2/3)*zp;
     plot3(xp,yp,zp,'+'); hold on;
  end
end

%  for m=3, put a point at the center of the slanted face

if(m==3)
 plot3(1/3,1/3,1/3,'o'); hold on; 
end

%---
% vertex coordinates
%----

x(1) = 0.0; y(1)= 0.0; z(1)= 0.0;
x(2) = 1.0; y(2)= 0.0; z(2)= 0.0;
x(3) = 0.0; y(3)= 1.0; z(3)= 0.0;
x(4) = 0.0; y(4)= 0.0; z(4)= 1.0;
                                                                                
% edge node labels

s(1,1)=2; s(1,2)=3; s(1,3)=4;
s(2,1)=3; s(2,2)=1; s(2,3)=4;
s(3,1)=1; s(3,2)=2; s(3,3)=4;
s(4,1)=1; s(4,2)=2; s(4,3)=3;

% draw the edges

for i=1:4
 for j=1:3
  xp(j)=x(s(i,j)); yp(j)=y(s(i,j)); zp(j)=z(s(i,j));

  % regular terahedron:
%     xp(j)=xp(j)+0.5*(yp(j)+zp(j)); yp(j)=0.5*sqrt(3)*(yp(j)+zp(j)/3.0);
%     zp(j) = sqrt(2/3)*zp(j);

 end
  plot3(xp,yp,zp,'k')
end

% formula = (m+1)*(m+2)*(m+3)/6.0

%---
% done
%---

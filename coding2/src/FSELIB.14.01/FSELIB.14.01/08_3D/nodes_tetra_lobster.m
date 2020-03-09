clear all
close all

%================================================
% Code nodes_tetra_lob
%
% Lobatto node distribution in the tetrahedron
% corresponding to a complete mth order expansion
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

figure(2)
hold on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)
set(gca,'fontsize',14)
axis('equal')
axis([0 1 0 1 0 1])
box on
view(56,22)

%----
% parameters
%----

m = 5

%-----
% draw the edges
%-----

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
  xpp(j)=x(s(i,j));
  ypp(j)=y(s(i,j));
  zpp(j)=z(s(i,j));
  % regular terahedron:
  xpr(j) =xpp(j)+0.5*(ypp(j)+zpp(j));
  ypr(j) =0.5*sqrt(3)*(ypp(j)+zpp(j)/3.0);
  zpr(j) = sqrt(2/3)*zpp(j);
 end
 figure(1)
 plot3(xpp,ypp,zpp,'k','linewidth',2)
 figure(2)
 plot3(xpr,ypr,zpr,'k','linewidth',2)
end

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
%---------------------------

Ic = 0; % node counter

%---
% xi-eta face
%---

for i=1:m+1;
 for j=1:m+2-i;
    Ic = Ic+1;
    l = m+3-i-j;
    xp = (1.0+2.0*v(i)-    v(j)-v(l))/3.0;
    yp = (1.0-    v(i)+2.0*v(j)-v(l))/3.0;
    zp = 0.0;
    figure(1)
    plot3(xp,yp,zp,'+');
    % regular tetrahedron:
    xp = xp+0.5*(yp+zp);
    yp= 0.5*sqrt(3)*(yp+zp/3.0);
    zp = sqrt(2/3)*zp;
    figure(2)
    plot3(xp,yp,zp,'+');
 end
end

%---
% eta-zeta face
%---

for j=1:m;
 for k=2:m+2-j;
    Ic = Ic+1;
    l = m+3-j-k;
    xp = 0.0;
    yp = (1.0+2.0*v(j)-    v(k)-v(l))/3.0;
    zp = (1.0-    v(j)+2.0*v(k)-v(l))/3.0;
    figure(1)
    plot3(xp,yp,zp,'o');
    % regular tetrahedron:
    xp = xp+0.5*(yp+zp);
    yp = 0.5*sqrt(3)*(yp+zp/3.0);
    zp = sqrt(2/3)*zp;
    figure(2)
    plot3(xp,yp,zp,'o');
 end
end

%---
% zeta-xi face
%---

for i=2:m;
 for k=2:m+2-i;
    Ic = Ic+1;
    l = m+3-i-k;
    xp = (1.0+2.0*v(i)-    v(k)-v(l))/3.0;
    yp = 0.0;
    zp = (1.0-    v(i)+2.0*v(k)-v(l))/3.0;
    figure(1)
    plot3(xp,yp,zp,'x');
    % regular tetrahedron:
    figure(2)
    xp=xp+0.5*(yp+zp);
    yp=0.5*sqrt(3)*(yp+zp/3.0);
    zp = sqrt(2/3)*zp;
    plot3(xp,yp,zp,'x');
 end
end

%---
% slanted face
%---

for i=2:m;
 for j=2:m+1-i;
    Ic = Ic+1;
    l = m+3-i-j;
    xp = (1.0+2.0*v(i)-     v(j)-v(l))/3.0;
    yp = (1.0-    v(i)+ 2.0*v(j)-v(l))/3.0;
    zp = 1.0-xp-yp;
    figure(1)
    plot3(xp,yp,zp,'p');
    % regular tetrahedron:
    figure(2)
    xp = xp+0.5*(yp+zp);
    yp = 0.5*sqrt(3)*(yp+zp/3.0);
    zp = sqrt(2/3)*zp;
    plot3(xp,yp,zp,'p');
 end
end

%---
% interior
%---

count = 0;

for i=2:m;
 for j=2:m+1-i
  for k=2:m+2-i-j

    l = m+4-i-j-k;
    xp = (1.0+3.0*v(i)-    v(j)-    v(k)-v(l))/4.0;
    yp = (1.0-    v(i)+3.0*v(j)-    v(k)-v(l))/4.0;
    zp = (1.0-    v(i)-    v(j)+3.0*v(k)-v(l))/4.0;
    figure(1)
    plot3(xp,yp,zp,'p','markersize',10,'color','red');
    % regular tetrahedron:
    xp = xp+0.5*(yp+zp);
    yp = 0.5*sqrt(3)*(yp+zp/3.0);
    zp = sqrt(2/3)*zp;
    figure(2)
    plot3(xp,yp,zp,'p','markersize',10,'color','red');
    Ic = Ic+1;
  end
 end
end

%---
% number of nodes
%---

Ic
formula = (m+1)*(m+2)*(m+3)/6.0



clear all
close all

%=======
% fekete_nodes
%
% prepare a graph of the element interpolation functions
% for the fekete node distribution
%=======

%---
% parameters
%---

m=6;

%---
% plotting
%---

figure(1)
hold on
axis equal
axis([0 1 0 1 0 1]);
axis off
set(gca,'fontsize',14)
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
box on

xx(1)=0.0;   yy(1)=0.0;
xx(2)=1.0;   yy(2)=0.0;
xx(3)=0.0;   yy(3)=1.0;
xx(4)=0.0;   yy(4)=0.0;

%---
% xihat etahat
%---

plot(xx,yy,'k');

figure(2)
hold on
axis equal
axis([0 1 0 1 0 1]);
axis off
set(gca,'fontsize',14)
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
box on

xx(1)=0.0;   yy(1)=0.0;
xx(2)=1.0;   yy(2)=0.0;
xx(3)=0.5;   yy(3)=sqrt(3)/2.0;
xx(4)=0.0;   yy(4)=0.0;

plot(xx,yy,'k')

%-------------
% lobatto grid
%-------------

[t, W] = lobatto(m-1);    % W is unused

v(1) = 0.0;
for i=2:m
 v(i) = 0.5*(1.0+t(i-1));
end
v(m+1) = 1.0;

for i=1:m+1
   for j=1:m+2-i
%     xx = v(i);
%     yy = v(j);
       k = m+3-i-j;
       xx = (1.0+2.0*v(i)-v(j)-v(k))/3.0;
       yy = (1.0+2.0*v(j)-v(i)-v(k))/3.0;
       figure(1)
        plot(xx,yy,'ko');
       figure(2)
        xx = xx+0.5*yy;
        yy = sqrt(3)*yy/2.0;
        plot(xx,yy,'ko')
   end
end

%-------------
% uniform grid
%-------------

for i=1:m+1
  v(i) = (i-1.0)/m;
end

for i=1:m+1
   for j=1:m+2-i
       xx = v(i);
       yy = v(j);
       figure(1)
        plot(xx,yy,'kx');
       figure(2)
        xx = xx+0.5*yy;
        yy = sqrt(3)*yy/2.0;
        plot(xx,yy,'kx')
   end
end

%=============
% fekete nodes
%=============

[nodes,xi,et,xih,eth] = fekete(m);

figure(1)
plot(xi,et,'k*','MarkerSize',10)

figure(2)
plot(xih,eth,'k*','MarkerSize',10)

%-----
% done
%-----

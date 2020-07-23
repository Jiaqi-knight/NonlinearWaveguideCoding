clear all
close all

%=====================
% triangulation driver
%=====================

ndiv = 0;
ndiv = 1;
ndiv = 2;
ndiv = 3;

a = 0.5;

%---
% triangulate
%---

[ne,ng,p,c,efl,gfl] = trgl6_sc(a, ndiv);

%---
% plot
%---

 figure(1)
 hold on
 axis equal
 xlabel('x','fontsize',15);
 ylabel('y','fontsize',15);
 set(gca,'fontsize',15)
 box on

 for l=1:ne
  j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2); 
  j=c(l,4); xp(2)=p(j,1); yp(2)=p(j,2); 
  j=c(l,2); xp(3)=p(j,1); yp(3)=p(j,2); 
  j=c(l,5); xp(4)=p(j,1); yp(4)=p(j,2); 
  j=c(l,3); xp(5)=p(j,1); yp(5)=p(j,2); 
  j=c(l,6); xp(6)=p(j,1); yp(6)=p(j,2); 
  j=c(l,1); xp(7)=p(j,1); yp(7)=p(j,2); 
  plot(xp, yp,'-ko'); 
 end

clear all
close all

%=================================
% triangular of an L-shaped domain
%=================================

ndiv = 3;
ndiv = 2;
ndiv = 0;
ndiv = 1;

%---
% triangulate
%---

[ne,ng,p,c,efl,gfl] = trgl6_L(ndiv);

%---
% plot
%---

 figure(1)
 hold on
 axis equal
 xlabel('x','fontsize',14);
 ylabel('y','fontsize',14);
 set(gca,'fontsize',14)
 box on

 for l=1:ne
  j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2); 
  j=c(l,4); xp(2)=p(j,1); yp(2)=p(j,2); 
  j=c(l,2); xp(3)=p(j,1); yp(3)=p(j,2); 
  j=c(l,5); xp(4)=p(j,1); yp(4)=p(j,2); 
  j=c(l,3); xp(5)=p(j,1); yp(5)=p(j,2); 
  j=c(l,6); xp(6)=p(j,1); yp(6)=p(j,2); 
  j=c(l,1); xp(7)=p(j,1); yp(7)=p(j,2); 
  plot(xp, yp,'-k'); 
 end

 for i=1:ng
  plot(p(i,1), p(i,2),'ko'); 
  if(gfl(i)==0)
%  plot(p(i,1), p(i,2),'ro'); 
  else
%  plot(p(i,1), p(i,2),'g*'); 
  end
 end

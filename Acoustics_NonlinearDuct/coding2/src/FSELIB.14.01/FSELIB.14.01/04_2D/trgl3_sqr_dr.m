clear all
close all

%=========================
% driver for triangulation
%=========================

ndiv = 0;
ndiv = 3;

%---
% triangulate
%---

[ne,ng,p,c,efl,gfl] = trgl3_sqr(ndiv);

%-----------------------------------------
% deform the square into a frivolous shape
%-----------------------------------------

ideform = 0;
ideform = 1;

if(ideform==1)

 for i=1:ng
  rad = sqrt(p(i,1)^2+p(i,2)^2);
  p(i,1) = p(i,1)*(0.25+0.50*rad);
  p(i,2) = p(i,2)*(1.00-0.10*rad);
 end

end

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
  j=c(l,2); xp(2)=p(j,1); yp(2)=p(j,2); 
  j=c(l,3); xp(3)=p(j,1); yp(3)=p(j,2); 
  j=c(l,1); xp(4)=p(j,1); yp(4)=p(j,2); 
  plot(xp, yp,'-ko'); 
 end


clear all
close all

%=============================================
% psi_tm
%
% prepare a graph of the element interpolation
% functions for a uniform grid on a triangle
%=============================================

figure(1)
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
set(gca,'fontsize',14)
box on

%----
% input parameters
%----

m = 4;
m = 19;

%---
% uniform master grid
%---

for i=1:m+1
 v(i) = (i-1.0)/m;
end


%---------------
% mark the nodes
%---------------

zz = 0.0

for i=1:m+1
  xx = v(i);
   for j=1:m+2-i
     yy = v(j);
     plot(xx,yy,'ko');
   end
end

nt = floor((m-1)/3)+1;

for q=1:nt
 x1 = (q-1)/m;
 y1 = x1;
 x2 = 1-(q-1)/m-x1;
 y2 = y1;
 x3 = x1;
 y3 = x2;
 x4 = x1;
 y4 = y1;
 plot([x1,x2,x3,x4],[y1,y2,y3,y4],'k-')
end

%-----
% done
%-----

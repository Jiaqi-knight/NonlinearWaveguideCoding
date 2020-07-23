clear all
close all

%======================================
% tetra4_sphere8
%
% driver for discretizing a sphere into
% 4-node tetrahedra
%======================================

%===
% prepare
%===

   figure(1)
   hold on
   axis equal
   xlabel('x','fontsize',14)
   ylabel('y','fontsize',14)
   zlabel('z','fontsize',14)
   set(gca,'fontsize',14)
   view([-26 14])
   box on
   axis([-1.25 1.25 -1.25 1.25 -1.25 1.25])
   axis off

%---
% parameters
%---

ndiv=0;
ndiv=2;
ndiv=3;
ndiv=1;

%---
% discretize
%---

[ne,ng,p,c,efl,gfl] = tetra4_sphere8 (ndiv);

%-------------------------
% randomize the interior nodes
%-------------------------

epsilon=0.10;
epsilon=0.00;

for i=1:ng
 if(gfl(i)==0)
  for j=1:3
   p(i,j) = p(i,j)+epsilon*(rand-0.5);
  end
 end
end

%---------
% plotting
%---------

option = 3;   % graded
option = 1;   % painted in color
option = 2;   % wire-frame 
option = 5;   % painted solid white

volume = 0.0;

for i=1:ne
  xmean = 0.0;
  ymean = 0.0;
  zmean = 0.0;
  for j=1:4
   xvis(j)=p(c(i,j),1);
   yvis(j)=p(c(i,j),2);
   zvis(j)=p(c(i,j),3);
   xmean = xmean+xvis(j);
   ymean = ymean+yvis(j);
   zmean = zmean+zvis(j);
  end
% if(xmean>0 | ymean>0)
   vol(i) = tetra4_vis (xvis,yvis,zvis,zvis,option);
% end
  volume = volume+vol(i);
end

volume

%---
% show the nodes
%---

for i=1:ng
  if(gfl(i)==1)
   plot3(p(i,1),p(i,2),p(i,3),'c.','MarkerSize',20)
  else
   plot3(p(i,1),p(i,2),p(i,3),'r.','MarkerSize',20)
  end
end

%-----
% done
%-----

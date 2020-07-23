%======================================
% tetra10_sphere8
%
% driver for discretizing a sphere into
% 10-node tetrahedra
%======================================

clear all
close all

%---
% parameters
%---

ndiv=3;
ndiv=0;
ndiv=1;
ndiv=2;

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
% discretize
%---

[ne,ng,p,c,efl,gfl] = tetra10_sphere8 (ndiv);

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

%---
% show the elements
%---

option = 2;   % wire-frame 
option = 1;   % painted in color
option = 3;   % graded
option = 5;   % painted solid white

volume = 0.0;

for i=1:ne
  xmean = 0.0;ymean = 0.0;zmean = 0.0;
  
  for j=1:4
   t = c(i,j);
   xvis(j)=p(t,1);
   yvis(j)=p(t,2);
   zvis(j)=p(t,3);
   xmean = xmean+xvis(j);
   ymean = ymean+yvis(j);
   zmean = zmean+zvis(j);
  end

  for j=5:10
   t = c(i,j);
   xvis(j)=p(t,1);
   yvis(j)=p(t,2);
   zvis(j)=p(t,3);
  end

%  if(xmean>0 | ymean>0)
   if(xmean>0)
   vol(i) = tetra10_vis (xvis,yvis,zvis,zvis,option);
   volume = volume+vol(i);
   end

  pause(0.00)

end

volume

%---
% plot the nodes
%---

inod = 0;

if(inod==1)

for i=1:ne
  for j=1:10
   t = c(i,j);
   if(gfl(t)==1)
    plot3(p(t,1),p(t,2),p(t,3),'c.','MarkerSize',20)
   else
    plot3(p(t,1),p(t,2),p(t,3),'r.','MarkerSize',20)
   end
  end
end

end

%-----
% done
%-----

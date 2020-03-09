function [ne,ng,p,c,efl,gfl] = trgl3_delaunay

%======================================
% Discretization of an arbitrary domain
% into 3-node triangular elements
% by the Delaunay triangulation
%======================================

%----------------------------------------
% Read the global nodes
% and boundary flags from file "points.dat"
%
% The three columns in this file are the
% x and y coordinates and the boundary flag
% of the global nodes:
%
%  x(1) y(1) gfl(1)
%  x(2) y(2) gfl(2)
%  ...  ...  ...
%  x(n) y(n) gfl(n)
%----------------------------------------

file1 = fopen('points.dat');

 points = fscanf(file1,'%f');

fclose(file1);

%---------------------------------------
% number, coordinates, and boundary flag
% of global nodes
%---------------------------------------

sp = size(points); ng=sp(1)/3;

for i=1:ng
 p(i,1) = points(3*i-2);
 p(i,2) = points(3*i-1);
 gfl(i) = points(3*i);
 xdel(i) = p(i,1);
 ydel(i) = p(i,2);
end

%-----------------------
% delaunay triangulation
%-----------------------

c = delaunay (xdel,ydel);

%-------------------------------
% extract the number of elements
%-------------------------------

sc = size(c); ne = sc(1,1);

%------------------------------------
% set the element-node boundary flags
%------------------------------------

for i=1:ne
 efl(i,1) = gfl(c(i,1));
 efl(i,2) = gfl(c(i,2));
 efl(i,3) = gfl(c(i,3));
end

%-----
% plot
%-----

figure(88)
hold on
axis equal
xlabel('x','fontsize',15);
ylabel('y','fontsize',15);
set(gca,'fontsize',15)
axis([0.35 0.85 0.25 0.65])

for i=1:ne
  xp(1) = p(c(i,1),1); yp(1) = p(c(i,1),2);
  xp(2) = p(c(i,2),1); yp(2) = p(c(i,2),2);
  xp(3) = p(c(i,3),1); yp(3) = p(c(i,3),2);
  xp(4) = p(c(i,1),1); yp(4) = p(c(i,1),2);
  plot(xp, yp,'k');
  if(efl(i,1)==1) plot(xp(1), yp(1), 'ro'); end;
  if(efl(i,2)==1) plot(xp(2), yp(2), 'r+'); end;
  if(efl(i,3)==1) plot(xp(3), yp(3), 'rx'); end;
end

%---------------------
% Voronoi tessellation
% (optional)
%---------------------

%figure(89);
 voronoi(xdel,ydel,':k')
 axis equal
 xlabel('x','fontsize',15);
 ylabel('y','fontsize',15);
 set(gca,'fontsize',15)
 axis([0.35 0.85 0.25 0.65])

%-----
% done
%-----

return

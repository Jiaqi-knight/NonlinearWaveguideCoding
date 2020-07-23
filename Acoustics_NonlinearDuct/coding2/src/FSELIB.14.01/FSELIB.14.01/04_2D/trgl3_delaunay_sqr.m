function [ne,ng,p,c,efl,gfl] = trgl3_delaunay_sqr

%===================================
% Delaunay triangulation of a square
%===================================

%---------------------
% window and grid size
%---------------------

X1 = -1.0; X2 = 1.0;
Y1 = -1.0; Y2 = 1.0;
N = 8; M = 8;

%--------
% prepare
%--------

Dx = (X2-X1)/N;
Dy = (Y2-Y1)/M;

%------------------------------
% arrange points on a mesh grid
% and set the boundary flag
%------------------------------

ng = 0;

for j=1:M+1
 for i=1:N+1

  ng = ng+1;
  p(ng,1) = X1+(i-1.0)*Dx;
  p(ng,2) = Y1+(j-1.0)*Dy;
  gfl(ng) = 0;
  if(i==1 | i==N+1 | j==1 | j==M+1) gfl(ng)=1; end;

 end
end

%-----------------------------
% randomize the interior nodes
%-----------------------------

Ic = N+2;

for j=2:M
 for i=2:N
   Ic=Ic+1;
   p(Ic,1) = p(Ic,1)+(rand-1.0)*0.5*Dx;
   p(Ic,2) = p(Ic,2)+(rand-1.0)*0.4*Dy;
 end
 Ic=Ic+2;
end

%-----------------------
% Delaunay triangulation
%-----------------------

for i=1:ng
 xdel(i) = p(i,1);
 ydel(i) = p(i,2);
end

c = delaunay(xdel,ydel);

%-------------------------------
% extract the number of elements
%-------------------------------

sc = size(c); ne = sc(1,1);

%-------------------------------
% set the element-node boundary flags
%-------------------------------

for i=1:ne
 efl(i,1) = gfl(c(i,1));
 efl(i,2) = gfl(c(i,2));
 efl(i,3) = gfl(c(i,3));
end

%--------------------------
% display the triangulation
%--------------------------

figure(88)
hold on
axis equal
xlabel('x','fontsize',15);
ylabel('y','fontsize',15);
set(gca,'fontsize',15)
axis([-1.1 1.1 -1.1 1.1])
box on

trimesh(c,p(:,1),p(:,2),zeros(ng,1),'EdgeColor','k');

%---------------------
% Voronoi tessellation
%---------------------

figure(88)
voronoi(xdel,ydel,':');
xlabel('x','fontsize',15);
ylabel('y','fontsize',15);
set(gca,'fontsize',15)

%----
% done
%-----

return

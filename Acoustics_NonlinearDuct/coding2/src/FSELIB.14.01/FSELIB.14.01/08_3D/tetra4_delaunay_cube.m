function [ne,ng,p,c,efl,gfl] = tetra4_delaunay_cube

%============================================
% Delaunay tessellation into 4-node tetrahedra
%============================================

%---------------------
% tessellation of a cube
% window and grid size
%---------------------

X1 = -1.0; X2 = 1.0;
Y1 = -1.0; Y2 = 1.0;
Z1 = -1.0; Z2 = 1.0;

Nx = 2; Ny = 2; Nz = 3;

%--------
% prepare
%--------

Dx = (X2-X1)/Nx;
Dy = (Y2-Y1)/Ny;
Dz = (Z2-Z1)/Nz;

%------------------------------
% arrange points on a mesh grid
% and set the boundary flag
%------------------------------

ng = 0;

for k=1:Nz+1
 for j=1:Ny+1
  for i=1:Nx+1
   ng = ng+1;
   p(ng,1) = X1+(i-1.0)*Dx;
   p(ng,2) = Y1+(j-1.0)*Dy;
   p(ng,3) = Z1+(k-1.0)*Dz;
   gfl(ng) = 0;
   if(i==1 | i==Nx+1 | j==1 | j==Ny+1 | k==1 | k==Nz+1) gfl(ng)=1; end;
   if(gfl(ng)==1)
     p(ng,1) = p(ng,1)+(rand-1.0)*0.00*Dx;
     p(ng,2) = p(ng,2)+(rand-1.0)*0.00*Dy;
     p(ng,3) = p(ng,3)+(rand-1.0)*0.00*Dz;
   end
  end
 end
end

%-----------------------
% Delaunay triangulation
%-----------------------

for i=1:ng
 xdel(i) = p(i,1);
 ydel(i) = p(i,2);
 zdel(i) = p(i,3);
end

c = delaunay3 (xdel,ydel,zdel);

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
 efl(i,4) = gfl(c(i,4));
end

%===
% ensure volumes are positive
%===

for i=1:ne

  for j=1:4
   xtmp(j)=p(c(i,j),1);
   ytmp(j)=p(c(i,j),2);
   ztmp(j)=p(c(i,j),3);
  end

  matrix=[1 1 1 1;
        xtmp(1) xtmp(2) xtmp(3) xtmp(4); ...
        ytmp(1) ytmp(2) ytmp(3) ytmp(4); ...
        ztmp(1) ztmp(2) ztmp(3) ztmp(4)];

  vlm = det(matrix)/6.0;

  if(vlm<0)
   save = c(i,1);
   c(i,1) = c(i,2);
   c(i,2) = save;
  end

end

%--------------------------
% display the triangulation
%--------------------------

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

%---------

option = 2;   % wire-frame 
option = 5;   % painted solid white
option = 3;   % graded
option = 1;   % painted in color

volume = 0.0;

for i=1:ne
  vol(i) = 0.0;
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

  if(xmean>0 | ymean>0)
% if(xmean>0)
   vol(i) = tetra4_vis (xvis,yvis,zvis,zvis,option);
   end
  volume = volume+vol(i);

end

volume

%===
% display the nodes
%===

for i=1:ng
  if(gfl(i) == 0)
%   plot3(p(i,1),p(i,2),p(i,3),'c.')
  elseif(gfl(i) == 1)
%   plot3(p(i,1),p(i,2),p(i,3),'r.')
  end
end

%----
% done
%-----

return

clear all
close all

%==============================
% tetra4_sub8_dr
%
% driver for subdividing
% a tetrahedron into 8 elements
%==============================

% constants
                                                                                
wait=0.5;

iopt = 1;
iopt = 2;

%---
if(iopt==1) % canonical vertex coordinates:
%---

  x(1) = 0.0; y(1) = 0.0; z(1)= 0.0;
  x(2) = 1.0; y(2) = 0.0; z(2)= 0.0;
  x(3) = 0.0; y(3) = 1.0; z(3)= 0.0;
  x(4) = 0.0; y(4) = 0.0; z(4)= 1.0;

 % canonical regular tetrahedron:

  for j=1:4
     x(j) = x(j)+0.5*(y(j)+z(j));
     y(j) = 0.5*sqrt(3)*(y(j)+z(j)/3.0);
     z(j) = sqrt(2/3)*z(j);
  end

%---
elseif(iopt==2) % arbitrary vertex coordinates:
%---
                                                                                
 x(1) = 0.0; y(1)= 0.1; z(1)= 0.0;
 x(2) = 0.9; y(2)= 0.2; z(2)= 0.5;
 x(3) = 0.1; y(3)= 1.1; z(3)=-0.05;
 x(4) =-0.1; y(4)=-0.1; z(4)= 1.2;

end

%---
% subdivide
%---

 [x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4 ...
 ,x5,y5,z5, x6,y6,z6, x7,y7,z7, x8,y8,z8 ] ...
 ...
 = tetra4_sub8(x,y,z);

%----------
% visualize
%----------

 f=x; f1=x; f2=x; f3=x; f4=x; f5=x; f6=x; f7=x; f8=x;

 ido=0;
 ido=1;

%---
 if(ido==1)
%---

 figure(88)
 hold on
 axis equal
 xlabel('x','fontsize',14)
 ylabel('y','fontsize',14)
 zlabel('z','fontsize',14)
 set(gca,'fontsize',14)
 view([55 34])
 box on
 axis off

 % display the parental element

 option = 2; % wireframe
 volume = tetra4_vis (x,y,z,f,option);

 % display the eight sub-elements

 option = 1; % in color

 for i=1:8
   vol(i) = 0.0;
 end

 vol(8) = tetra4_vis (x8, y8, z8, f8, option);pause(wait);
 vol(7) = tetra4_vis (x7, y7, z7, f7, option);pause(wait);
 vol(6) = tetra4_vis (x6, y6, z6, f6, option);pause(wait);
 vol(5) = tetra4_vis (x5, y5, z5, f5, option);pause(wait);
 vol(4) = tetra4_vis (x4, y4, z4, f4, option);pause(wait);
 vol(3) = tetra4_vis (x3, y3, z3, f3, option);pause(wait);
 vol(2) = tetra4_vis (x2, y2, z2, f2, option);pause(wait);
 vol(1) = tetra4_vis (x1, y1, z1, f1, option);pause(wait);

%---
 end
%---

%--------
% volumes
%--------

 vlm= 0.0;
 for i=1:8
%   vol(i)
  vlm = vlm+vol(i);
 end

[vlm volume] % should be the same

%---------------------
% display the individual subelements
%---------------------

 ido=0;
 ido=1;

%---
 if(ido==1)
%---

 for i=1:8
   figure(i)
   hold on
   axis equal
   xlabel('x','fontsize',14)
   ylabel('y','fontsize',14)
   zlabel('z','fontsize',14)
   set(gca,'fontsize',14)
   view([56 16])
   box on
%  axis([0 1 0 1 0 1])
   axis off
   option = 2; % wireframe
   volume = tetra4_vis (x,y,z,f,option);
   option = 1; % in color
   if(i==1)
     vol(1) = tetra4_vis (x1, y1, z1, f1, option);
   elseif(i==2)
     vol(2) = tetra4_vis (x2, y2, z2, f2, option);
   elseif(i==3)
     vol(3) = tetra4_vis (x3, y3, z3, f3, option);
   elseif(i==4)
     vol(4) = tetra4_vis (x4, y4, z4, f4, option);
   elseif(i==5)
     vol(5) = tetra4_vis (x5, y5, z5, f5, option);
   elseif(i==6)
     vol(6) = tetra4_vis (x6, y6, z6, f6, option);
   elseif(i==7)
     vol(7) = tetra4_vis (x7, y7, z7, f7, option);
   elseif(i==8)
     vol(8) = tetra4_vis (x8, y8, z8, f8, option);
   end
 end

%---
 end
%---

%-----
% done
%-----

%========================
% tetra4
%
% visualize a tetrahedron
%========================

close all
clear all

%===
% prepare
%===

figure(1)
hold on;
xlabel('\Xi','fontsize',14)
ylabel('H','fontsize',14)
zlabel('Z','fontsize',14)
set(gca,'fontsize',14)
axis('equal')
box on
axis([-0.3 1.3 -0.3 1.3 -0.3 1.3]);
view([62 10])

%---
imenu = 1;
imenu = 2;
%---

% face node labels:

 s(1,1)=3; s(1,2)=2; s(1,3)=4;    % face 1
 s(2,1)=1; s(2,2)=3; s(2,3)=4;    % face 2
 s(3,1)=1; s(3,2)=4; s(3,3)=2;    % face 3
 s(4,1)=1; s(4,2)=2; s(4,3)=3;    % face 4

% arbitrary vertex coordinates:

 if(imenu==1)
  x(1) = 0.0; y(1)= 0.0; z(1)= 0.0;
  x(2) = 0.9; y(2)= 0.2; z(2)= 0.5;
  x(3) = 0.1; y(3)= 1.1; z(3)=-0.05;
  x(4) =-0.1; y(4)=-0.1; z(4)= 1.2;
 end

% canonical vertex coordinates:

 if(imenu==2)
  x(1) = 0.0; y(1) = 0.0; z(1)= 0.0;
  x(2) = 1.0; y(2) = 0.0; z(2)= 0.0;
  x(3) = 0.0; y(3) = 1.0; z(3)= 0.0;
  x(4) = 0.0; y(4) = 0.0; z(4)= 1.0;
 end

%---
% paint the faces:
%---

option=2; % wireframe
option=1; % solid

for i=1:4

 for j=1:3

  xv = x(s(i,j));
  yv = y(s(i,j));
  zv = z(s(i,j));

  xp(j)=xv; yp(j)=yv; zp(j)=zv;

  if(imenu==2)%  canonical regular tetrahedron:
   xp(j)=xv+0.5*(yv+zv);
   yp(j)=0.5*sqrt(3)*(yv+zv/3.0); 
   zp(j) = sqrt(2/3)*zv;
  end

 end
    if(option==1) % solid
       if(i==1) patch(xp,yp,zp,'r'); end
       if(i==2) patch(xp,yp,zp,'b'); end
       if(i==3) patch(xp,yp,zp,'y'); end
       if(i==4) patch(xp,yp,zp,'g'); end
     else    % wireframe
       plot3 (xp,yp,zp); hold on;
     end
end
                                                                                
%----------------------------------------------
% gradients of the node interpolation functions
%----------------------------------------------

x21=x(2)-x(1); y21=y(2)-y(1); z21=z(2)-z(1);
x31=x(3)-x(1); y31=y(3)-y(1); z31=z(3)-z(1);
x41=x(4)-x(1); y41=y(4)-y(1); z41=z(4)-z(1);
x32=x(3)-x(2); y32=y(3)-y(2); z32=z(3)-z(2);
x42=x(4)-x(2); y42=y(4)-y(2); z42=z(4)-z(2);

%----------------
% Jacobian matrix
%----------------

Jac(1,1)=x21; Jac(1,2)=x31; Jac(1,3)=x41;
Jac(2,1)=y21; Jac(2,2)=y31; Jac(2,3)=y41;
Jac(3,1)=z21; Jac(3,2)=z31; Jac(3,3)=z41;

%-------
% volume
%-------

V = det(Jac)/6.0;

V6 = 6.0*V;

%---
% first gradient 
%---

gx(1) = (-y32*z42 + z32*y42)/V6;
gy(1) = ( x32*z42 - z32*x42)/V6;
gz(1) = (-x32*y42 + y32*x42)/V6;

mag = sqrt(gx(1)^2+gy(1)^2+gz(1)^2);
fc = -0.10/mag;
px(1) = (x(2)+x(3)+x(4))/3.0;
py(1) = (y(2)+y(3)+y(4))/3.0;
pz(1) = (z(2)+z(3)+z(4))/3.0;
%  quiver3(px(1),py(1),pz(1),fc*gx(1),fc*gy(1),fc*gz(1),'ob')
px(2)=px(1)+fc*gx(1);
py(2)=py(1)+fc*gy(1);
pz(2)=pz(1)+fc*gz(1);
plot3(px,py,pz,'.-k')

%---
% second gradient 
%---

gx(2) = ( y31*z41 - z31*y41)/V6;
gy(2) = (-x31*z41 + z31*x41)/V6;
gz(2) = ( x31*y41 - y31*x41)/V6;

mag = sqrt(gx(2)^2+gy(2)^2+gz(2)^2);
fc = -0.10/mag;
px(1) = (x(1)+x(3)+x(4))/3.0;
py(1) = (y(1)+y(3)+y(4))/3.0;
pz(1) = (z(1)+z(3)+z(4))/3.0;
%quiver3(px(1),py(1),pz(1),fc*gx(2),fc*gy(2),fc*gz(2),'or')
px(2)=px(1)+fc*gx(2);
py(2)=py(1)+fc*gy(2);
pz(2)=pz(1)+fc*gz(2);

%plot3(px,py,pz,'.-k')

%---
% third gradient 
%---

gx(3) = (-y21*z41 + z21*y41)/V6;
gy(3) = ( x21*z41 - z21*x41)/V6;
gz(3) = (-x21*y41 + y21*x41)/V6;

mag = sqrt(gx(3)^2+gy(3)^2+gz(3)^2);
fc = -0.10/mag;
px(1) = (x(1)+x(2)+x(4))/3.0;
py(1) = (y(1)+y(2)+y(4))/3.0;
pz(1) = (z(1)+z(2)+z(4))/3.0;
%quiver3(px(1),py(1),pz(1),fc*gx(3),fc*gy(3),fc*gz(3),'og')
px(2)=px(1)+fc*gx(3);
py(2)=py(1)+fc*gy(3);
pz(2)=pz(1)+fc*gz(3);

%plot3(px,py,pz,'.-k')

%---
% fourth gradient 
%---

gx(4) = ( y21*z31 - z21*y31)/V6;
gy(4) = (-x21*z31 + z21*x31)/V6;
gz(4) = ( x21*y31 - y21*x31)/V6;

mag = sqrt(gx(4)^2+gy(4)^2+gz(4)^2);
fc = -0.10/mag;
px(1) = (x(1)+x(2)+x(3))/3.0;
py(1) = (y(1)+y(2)+y(3))/3.0;
pz(1) = (z(1)+z(2)+z(3))/3.0;
%quiver3(px(1),py(1),pz(1),fc*gx(4),fc*gy(4),fc*gz(4),'ok')
px(2) = px(1)+fc*gx(4);
py(2) = py(1)+fc*gy(4);
pz(2) = pz(1)+fc*gz(4);

%plot3(px,py,pz,'.-k')

%-----
% done
%-----

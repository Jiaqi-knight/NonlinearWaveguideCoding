function volume = tetra4_vis (x,y,z,f,option)

%========================
% tetra4_vis
%
% visualize a tetrahedron defined by
% four vertices hosted in  x, y, z
%========================

%---
% face node labels
%---

s(1,1)=2; s(1,2)=3; s(1,3)=4;    % face 1
s(2,1)=3; s(2,2)=1; s(2,3)=4;    % face 2
s(3,1)=1; s(3,2)=2; s(3,3)=4;    % face 3
s(4,1)=1; s(4,2)=2; s(4,3)=3;    % face 4

%---
% paint the faces
%---

for i=1:4

 for j=1:3
  xp(j) = x(s(i,j));
  yp(j) = y(s(i,j));
  zp(j) = z(s(i,j));
  fp(j) = f(s(i,j));
 end

 xp(4)=xp(1); yp(4)=yp(1); zp(4)=zp(1); fp(4)=fp(1);

  if(option==1)   % solid

    if(i==1) patch (xp,yp,zp,'r'); end
    if(i==2) patch (xp,yp,zp,'b'); end
    if(i==3) patch (xp,yp,zp,'y'); end
    if(i==4) patch (xp,yp,zp,'g'); end

  elseif(option==2)    % wireframe

    plot3 (xp,yp,zp,'k');

  elseif(option==3)  % graded

    cp(1)=0.5+0.3*xp(1);
    cp(2)=0.5+0.3*xp(2);
    cp(3)=0.5+0.3*xp(3);
    cp(4)=0.5+0.3*xp(4);

    cp(1)=0.5+0.3*zp(1);
    cp(2)=0.5+0.3*zp(2);
    cp(3)=0.5+0.3*zp(3);
    cp(4)=0.5+0.3*zp(4);

    patch(xp,yp,zp,cp);

   elseif(option==4)  % graded

    patch(xp,yp,zp,fp);

   elseif(option==5)   % solid white

    patch(xp,yp,zp,'w');

   end
end

%-------------------
% compute the volume
%-------------------

matrix=[1 1 1 1;
        x(1) x(2) x(3) x(4); ...
        y(1) y(2) y(3) y(4); ...
        z(1) z(2) z(3) z(4)];

volume = det(matrix)/6.0;

%-----
% done
%-----

return

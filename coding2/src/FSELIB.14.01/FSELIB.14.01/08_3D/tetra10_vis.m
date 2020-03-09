function volume = tetra10_vis (x,y,z,f,option)

%========================
% tetra10_vis
%
% visualize a ten-node tetrahedron
%========================

%---
% paint the faces:
%---

for i=1:16
 
  if(i==1) 
    xp(1)=x(1); yp(1)=y(1); zp(1)=z(1); fp(1)=f(1);
    xp(2)=x(5); yp(2)=y(5); zp(2)=z(5); fp(2)=f(5);
    xp(3)=x(7); yp(3)=y(7); zp(3)=z(7); fp(3)=f(7);
    xp(4)=x(1); yp(4)=y(1); zp(4)=z(1); fp(4)=f(1);
  elseif(i==2) 
    xp(1)=x(2); yp(1)=y(2); zp(1)=z(2); fp(1)=f(2);
    xp(2)=x(6); yp(2)=y(6); zp(2)=z(6); fp(2)=f(6);
    xp(3)=x(5); yp(3)=y(5); zp(3)=z(5); fp(3)=f(5);
    xp(4)=x(2); yp(4)=y(2); zp(4)=z(2); fp(4)=f(2);
  elseif(i==3) 
    xp(1)=x(3); yp(1)=y(3); zp(1)=z(3); fp(1)=f(3);
    xp(2)=x(7); yp(2)=y(7); zp(2)=z(7); fp(2)=f(7);
    xp(3)=x(6); yp(3)=y(6); zp(3)=z(6); fp(3)=f(6);
    xp(4)=x(3); yp(4)=y(3); zp(4)=z(3); fp(4)=f(3);
  elseif(i==4) 
    xp(1)=x(5); yp(1)=y(5); zp(1)=z(5); fp(1)=f(5);
    xp(2)=x(6); yp(2)=y(6); zp(2)=z(6); fp(2)=f(6);
    xp(3)=x(7); yp(3)=y(7); zp(3)=z(7); fp(3)=f(7);
    xp(4)=x(5); yp(4)=y(5); zp(4)=z(5); fp(4)=f(5);
  elseif(i==5) 
    xp(1)=x(1); yp(1)=y(1); zp(1)=z(1); fp(1)=f(1);
    xp(2)=x(5); yp(2)=y(5); zp(2)=z(5); fp(2)=f(5);
    xp(3)=x(8); yp(3)=y(8); zp(3)=z(8); fp(3)=f(8);
    xp(4)=x(1); yp(4)=y(1); zp(4)=z(1); fp(4)=f(1);
  elseif(i==6) 
    xp(1)=x(2); yp(1)=y(2); zp(1)=z(2); fp(1)=f(2);
    xp(2)=x(9); yp(2)=y(9); zp(2)=z(9); fp(2)=f(9);
    xp(3)=x(5); yp(3)=y(5); zp(3)=z(5); fp(3)=f(5);
    xp(4)=x(2); yp(4)=y(2); zp(4)=z(2); fp(4)=f(2);
  elseif(i==7) 
    xp(1)=x(4); yp(1)=y(4); zp(1)=z(4); fp(1)=f(4);
    xp(2)=x(8); yp(2)=y(8); zp(2)=z(8); fp(2)=f(8);
    xp(3)=x(9); yp(3)=y(9); zp(3)=z(9); fp(3)=f(9);
    xp(4)=x(4); yp(4)=y(4); zp(4)=z(4); fp(4)=f(4);
  elseif(i==8) 
    xp(1)=x(5); yp(1)=y(5); zp(1)=z(5); fp(1)=f(5);
    xp(2)=x(9); yp(2)=y(9); zp(2)=z(9); fp(2)=f(9);
    xp(3)=x(8); yp(3)=y(8); zp(3)=z(8); fp(3)=f(8);
    xp(4)=x(5); yp(4)=y(5); zp(4)=z(5); fp(4)=f(5);
  elseif(i==9) 
    xp(1)=x(1); yp(1)=y(1); zp(1)=z(1); fp(1)=f(1);
    xp(2)=x(8); yp(2)=y(8); zp(2)=z(8); fp(2)=f(8);
    xp(3)=x(7); yp(3)=y(7); zp(3)=z(7); fp(3)=f(7);
    xp(4)=x(1); yp(4)=y(1); zp(4)=z(1); fp(4)=f(1);
  elseif(i==10) 
    xp(1)=x(8);  yp(1)=y(8);  zp(1)=z(8);  fp(1)=f(8);
    xp(2)=x(4);  yp(2)=y(4);  zp(2)=z(4);  fp(2)=f(4);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(8);  yp(4)=y(8);  zp(4)=z(8);  fp(4)=f(8);
  elseif(i==11) 
    xp(1)=x(3);  yp(1)=y(3);  zp(1)=z(3);  fp(1)=f(3);
    xp(2)=x(7);  yp(2)=y(7);  zp(2)=z(7);  fp(2)=f(7);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(3);  yp(4)=y(3);  zp(4)=z(3);  fp(4)=f(3);
  elseif(i==12) 
    xp(1)=x(7);  yp(1)=y(7);  zp(1)=z(7);  fp(1)=f(7);
    xp(2)=x(8);  yp(2)=y(8);  zp(2)=z(8);  fp(2)=f(8);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(7);  yp(4)=y(7);  zp(4)=z(7);  fp(4)=f(7);
  elseif(i==13) 
    xp(1)=x(2); yp(1)=y(2); zp(1)=z(2); fp(1)=f(2);
    xp(2)=x(6); yp(2)=y(6); zp(2)=z(6); fp(2)=f(6);
    xp(3)=x(9); yp(3)=y(9); zp(3)=z(9); fp(3)=f(9);
    xp(4)=x(2); yp(4)=y(2); zp(4)=z(2); fp(4)=f(2);
  elseif(i==14) 
    xp(1)=x(6);  yp(1)=y(6);  zp(1)=z(6);  fp(1)=f(6);
    xp(2)=x(3);  yp(2)=y(3);  zp(2)=z(3);  fp(2)=f(3);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(6);  yp(4)=y(6);  zp(4)=z(6);  fp(4)=f(6);
  elseif(i==15) 
    xp(1)=x(4);  yp(1)=y(4);  zp(1)=z(4);  fp(1)=f(4);
    xp(2)=x(9);  yp(2)=y(9);  zp(2)=z(9);  fp(2)=f(9);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(4);  yp(4)=y(4);  zp(4)=z(4);  fp(4)=f(4);
  elseif(i==16) 
    xp(1)=x(9);  yp(1)=y(9);  zp(1)=z(9);  fp(1)=f(9);
    xp(2)=x(6);  yp(2)=y(6);  zp(2)=z(6);  fp(2)=f(6);
    xp(3)=x(10); yp(3)=y(10); zp(3)=z(10); fp(3)=f(10);
    xp(4)=x(9);  yp(4)=y(9);  zp(4)=z(9);  fp(4)=f(9);
  end

  if(option==1)   % solid

    if(i==1)  patch(xp,yp,zp,'r'); end
    if(i==2)  patch(xp,yp,zp,'b'); end
    if(i==3)  patch(xp,yp,zp,'y'); end
    if(i==4)  patch(xp,yp,zp,'g'); end
    if(i==5)  patch(xp,yp,zp,'r'); end
    if(i==6)  patch(xp,yp,zp,'b'); end
    if(i==7)  patch(xp,yp,zp,'y'); end
    if(i==8)  patch(xp,yp,zp,'g'); end
    if(i==9)  patch(xp,yp,zp,'r'); end
    if(i==10) patch(xp,yp,zp,'b'); end
    if(i==11) patch(xp,yp,zp,'y'); end
    if(i==12) patch(xp,yp,zp,'g'); end
    if(i==13) patch(xp,yp,zp,'r'); end
    if(i==14) patch(xp,yp,zp,'b'); end
    if(i==15) patch(xp,yp,zp,'y'); end
    if(i==16) patch(xp,yp,zp,'g'); end

  elseif(option==2)    % wireframe

    plot3(xp,yp,zp,'k');

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

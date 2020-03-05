function [X,Y,Z] = Tubefxn(Rfxn,Nr,Paramtricfxn,tbound,Nt)
% This function is the numerial analogy to the tube function in maple. It
% generates the X,Y and Z coordinates to be used by the surf function to
% generate the intended tube figure

% Rfxn is the parametric function of Radius as a function of t
% Paramtricfxn is a cell of path parametric functions
% tbound range is the range of values of t
% Nr is the number of steps in r
% Nt is the number of steps in t

Example1:
[X,Y,Z] = Tubefxn(@(t) t/3,40,{@(t) t.*sin(t),@(t) t.*cos(t),@(t) t},[0,20],200);
s = surf(X,Y,Z); set(s,'FaceColor','r'); axis equal
set(gca,'Xdir','reverse','Zdir','reverse')
camlight 
lighting gouraud
title('tubeN(@(t) t/3,40,{@(t) t.*sin(t),@(t) t.*cos(t),@(t) t},[0,10],20)')

Example2:
[X,Y,Z] = Tubefxn(@(t) sin(t),20,{@(t) t,@(t) 0*t,@(t) 0*t},[0,10],40);
s = surf(X,Y,Z); set(s,'FaceColor','r'); axis equal
camlight 
lighting gouraud
title('Tubefxn(@(t) sin(t),20,{@(t) t,@(t) 0*t,@(t) 0*t},[0,10],40)')

Example3:
[X,Y,Z] = Tubefxn(@(t) t,40,{@(t) 0*t,@(t) 0*t,@(t) sin(t)},[0,10],20);
s = surf(X,Y,Z); set(s,'FaceColor','r'); axis equal
camlight 
lighting gouraud
title('Tubefxn(@(t) t,40,{@(t) 0*t,@(t) 0*t,@(t) sin(t)},[0,10],20)')
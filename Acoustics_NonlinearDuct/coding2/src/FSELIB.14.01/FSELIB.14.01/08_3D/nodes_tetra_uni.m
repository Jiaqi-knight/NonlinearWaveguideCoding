clear all
close all

%=================================================
% uniform node distribution in the tetrahedron
% corresponding to a complete mth order expansion
%================================================

%----
% prepare
%----

figure(1)
hold on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)
set(gca,'fontsize',14)
axis('equal')
axis([ 0 1 0 1 0 1])
box on
view(56,22)

%----
% parameters
%----

m = 5

%--------------------
% uniform master grid
%--------------------

for i=1:m+1
   v(i) = (i-1.0)/m;
end

%---------------------------
% deploy and count the nodes
%---------------------------

count = 0;

for i=1:m+1;
 xp = v(i);
 for j=1:m+2-i
 yp = v(j);
  for k=1:m+3-i-j
   zp = v(k);

    count = count+1;

    l = m+4-i-j-k;
    if(l==1) plot3(xp,yp,zp,'k+'); end;
    if(l==2) plot3(xp,yp,zp,'ko'); end;
    if(l==3) plot3(xp,yp,zp,'kx'); end;
    if(l==4) plot3(xp,yp,zp,'k*'); end;
    if(l==5) plot3(xp,yp,zp,'k+'); end;
    if(l==6) plot3(xp,yp,zp,'ko'); end;
    if(l==7) plot3(xp,yp,zp,'kx'); end;
    if(l==8) plot3(xp,yp,zp,'k*'); end;
    if(l==9) plot3(xp,yp,zp,'k+'); end;

  end
 end
end

x(1) = 0.0; y(1)= 0.0; z(1)= 0.0;
x(2) = 1.0; y(2)= 0.0; z(2)= 0.0;
x(3) = 0.0; y(3)= 1.0; z(3)= 0.0;
x(4) = 0.0; y(4)= 0.0; z(4)= 1.0;
                                                                                
% edge node labels

s(1,1)=2; s(1,2)=3; s(1,3)=4;
s(2,1)=3; s(2,2)=1; s(2,3)=4;
s(3,1)=1; s(3,2)=2; s(3,3)=4;
s(4,1)=1; s(4,2)=2; s(4,3)=3;
                                                                                
% draw the edges

for i=1:4
 for j=1:3
  xp(j)=x(s(i,j)); yp(j)=y(s(i,j)); zp(j)=z(s(i,j));
  plot3(xp,yp,zp,'k')
 end
end

count
formula = (m+1)*(m+2)*(m+3)/6.0


for i=2:m
 q = v(i);
 plot3([q,0,0,q],[0,q,0,0],[0,0,q,0],':k')
end

%-----
% done
%-----

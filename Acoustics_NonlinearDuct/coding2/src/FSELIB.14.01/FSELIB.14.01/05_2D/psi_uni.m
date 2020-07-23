
clear all
close all

%=============================================
% prepare a graph of the element interpolation
% functions for a uniform grid on a triangle
%=============================================

%---
% prepare
%---

figure(1)
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
box on
view([56,30])

%----
% input parameters
%----

m = 4;
m = 3;

in = 1; jn = 5;   % choose the node (i,j)
in = 3; jn = 2;   % choose the node (i,j)
in = 2; jn = 2;   % choose the node (i,j)

%---
% uniform master grid
%---

for i=1:m+1
 v(i) = (i-1.0)/m;
end

%---
% prepare
%---

kn = m+3-in-jn;

%---
% define the grid lines for plotting purposes
%---

N=32; M=32;

Dx=1.0/N; Dy=1.0/M;

for i=1:N+1 
  x(i)=Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end

for j=1:M+1
  y(j)=Dy*(j-1.0);
% fprintf (1,'%6.2f',Y(j))
end

%------------------------------
% loop over vertical grid lines
%------------------------------

for i=1:N
for j=1:M+2-i

 xp(j) = x(i);
 yp(j) = y(j);
 zeta= 1.0-xp(j)-yp(j);

%---
% compute the interpolation function: zp
%---

 zp(j) = 1.0;
 if(in>1) 
 for ic=1:in-1
  zp(j) = zp(j)*(xp(j)-v(ic))/(v(in)-v(ic));
 end
 end

 if(jn>1) 
 for jc=1:jn-1
  zp(j) = zp(j)*(yp(j)-v(jc))/(v(jn)-v(jc));
 end
 end

 if(kn>1) 
 for kc=1:kn-1
  zp(j) = zp(j)*(zeta-v(kc))/(v(kn)-v(kc));
 end
 end

end

plot3(xp,yp,zp,'k'); hold on;

end

clear xp,yp,zp

%--------------------------------
% loop over horizontal grid lines
%--------------------------------

for j=1:M
for i=1:N+2-j

 yp(i) = y(j);
 xp(i) = x(i);
 zeta= 1.0-xp(i)-yp(i);

%---
% compute the interpolation function: zp
%---

 zp(i) = 1.0;
 if(in>1) 
 for ic=1:in-1
  zp(i) = zp(i)*(xp(i)-v(ic))/(v(in)-v(ic));
 end
 end
 if(jn>1) 
 for jc=1:jn-1
  zp(i) = zp(i)*(yp(i)-v(jc))/(v(jn)-v(jc));
 end
 end
 if(kn>1) 
 for kc=1:kn-1
  zp(i) = zp(i)*(zeta-v(kc))/(v(kn)-v(kc));
 end
 end
 
end

plot3(xp,yp,zp,'k'); hold on;

end

%------------------------------------
% plot the hypotenese of the triangle
%------------------------------------

xx(1)=0.0;yy(1)=0.0;zz(1)=0.0;
xx(2)=1.0;yy(2)=0.0;zz(2)=0.0;
xx(3)=0.0;yy(3)=1.0;zz(3)=0.0;
xx(4)=0.0;yy(4)=0.0;zz(4)=0.0;

plot3(xx,yy,zz,'k-');

%---------------
% mark the nodes
%---------------

zz = 0.0

for i=1:m+1
  xx = v(i);
   for j=1:m+2-i
     yy = v(j);
     plot3 (xx,yy,zz,'ko');
   end
end

%-----
% done
%-----

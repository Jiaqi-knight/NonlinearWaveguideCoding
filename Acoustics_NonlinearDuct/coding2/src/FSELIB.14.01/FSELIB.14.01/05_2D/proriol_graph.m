
clear all
close all

%===========================================
% proriol_graph
%
% prepare graphs of the Proriol polynomials
%===========================================

%-----------
% parameters
%-----------

k = 2;  % order
l = 0;  % order

mode = 11;

%---
% prepare
%---

figure(1)
hold on
%axis equal
set(gca,'fontsize',14)
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('P','fontsize',14);
view(40,26)
if(k==2 & l==0)
 axis([0 1 0 1 -0.5 1])
else
 axis([0 1 0 1 -1 1])
end
box on

%--------------
% plotting grid
%--------------

N=32; M=32;   % plotting grid

Dx=1.0/N; Dy=1.0/M;

for i=1:N+1 
  x(i)=Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end

for j=1:M+1
  y(j)=Dy*(j-1.0);
% fprintf (1,'%6.2f',Y(j))
end

%---
% vertical lines
%---

for i=1:N
for j=1:M+2-i
 xp(j) = x(i);
 yp(j) = y(j);
 zeta= 1.0-xp(j)-yp(j);

 xii = xp(j); eta=yp(j);
 if(eta>0.999999) eta=0.999999; end
 xip = 2*xii/(1-eta)-1; etap = 2*eta-1;

 if(mode==10) 
   pr = xii-zeta; % P_10
 elseif(mode==01) 
   pr = 3*eta-1;   % P_01
 elseif(mode==20) 
   pr = 6*xii^2 + 6*xii*eta + eta^2 - 6*xii-2*eta + 1; % P_20
 elseif(mode==11) 
   pr = (2*xii+eta-1)*(5*eta-1);  % P_11
 elseif(mode==02) 
   pr = 10*eta^2 - 8*eta + 1; % P_02
 end

 zp(j) = pr;
 zp(j) = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;
% zp(j) = zp(j)^2;

end
plot3(xp,yp,zp,'k')
end

clear xp yp zp

%---
% horizontal lines
%---

for j=1:M
for i=1:N+2-j
 yp(i) = y(j);
 xp(i) = x(i);
 zeta= 1.0-xp(i)-yp(i);

 xii = xp(i); eta=yp(i);
 if(eta>0.999999) eta=0.999999; end
 xip = 2*xii/(1-eta)-1; etap = 2*eta-1;

 if(mode==10) 
   pr = xii-zeta; % P_10
 elseif(mode==01) 
   pr = 3*eta-1;   % P_01
 elseif(mode==20) 
   pr = 6*xii^2 + 6*xii*eta + eta^2 - 6*xii-2*eta + 1; % P_20
 elseif(mode==11) 
   pr = (2*xii+eta-1)*(5*eta-1);  % P_11
 elseif(mode==02) 
   pr = 10*eta^2 - 8*eta + 1; % P_02
 end

 zp(i) = pr;
 zp(i) = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;

% zp(i) = zp(i)^2;

end
plot3(xp,yp,zp,'k');
end

%------------------
% plot the triangle
%------------------

xx(1)=0.0; yy(1)=0.0; zz(1)=0.0;
xx(2)=1.0; yy(2)=0.0; zz(2)=0.0;
xx(3)=0.0; yy(3)=1.0; zz(3)=0.0;
xx(4)=0.0; yy(4)=0.0; zz(4)=0.0;

plot3(xx,yy,zz,'k');

%-----
% done
%-----

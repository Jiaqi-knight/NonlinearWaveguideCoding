
clear all
close all

%===========================================
% psi_lob
%
% graphs of the interpolation functions
% for the lobatto triangle nodes
%===========================================

%---
% parameters
%---

m = 7;  % polynomial order

%---
% choose the node
% nodes are counted in vertical layers
%---

node =12;  % node
node =15;  % node
node = 1;  % node
node =22;  % node

%---
% prepare
%---

figure(1)
hold on
% axis equal
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
view([55,54])
box on

%-------------
% lobatto grid
%-------------

[t, W] = lobatto(m-1);    % W will not be needed

v(1) = 0.0;

for i=2:m
 v(i) = 0.5*(1.0+t(i-1));
end

v(m+1) = 1.0;

%-------------------------------------------
% compute the generalized vandermonde matrix
% based on the proriol polynomials
% nodes counted with in vertical layers
% by the indexes Ic and Jc
%-------------------------------------------

Ic = 0;

for i=1:m+1
 for j=1:m+2-i

  Ic = Ic+1;

  k = m+3-i-j;
  xi  = (1.0+ 2*v(i)  -v(j)-v(k))/3.0;
  eta = (1.0-   v(i)+2*v(j)-v(k))/3.0;
    if(eta>0.999999) eta=0.999999; end;

  xip = 2*xi/(1-eta)-1;
  etap = 2*eta-1;

  Jc=0;
  for k=0:m
   for l=0:m-k
    Jc=Jc+1;
    van(Jc,Ic) = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;
   end
  end

  plot3(xi,eta,0,'ko');

  end
end

%---
% inverse of the Vandermond matrix
%---

vani = inv(van);

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

%---------------
% vertical lines
%---------------

for i=1:N
 for j=1:M+2-i

 xp(j) = x(i);
 yp(j) = y(j);
 zeta= 1.0-xp(j)-yp(j);

%---
 xii = xp(j); eta=yp(j);
 if(eta>0.999999) eta=0.999999; end
 xip = 2*xii/(1-eta)-1; etap = 2*eta-1;
%---

  Jc=0;
  for k=0:m
   for l=0:m-k
    Jc=Jc+1;
    pror(Jc) = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;
   end
  end

  psil = vani*pror';

  zp(j) = psil(node);

 end
plot3(xp,yp,zp,'k')
end

clear xp yp zp

%-----------------
% horizontal lines
%-----------------

for j=1:M
for i=1:N+2-j
 yp(i) = y(j);
 xp(i) = x(i);

%---
 xii = xp(i);
 eta = yp(i);
 if(eta>0.999999) eta=0.999999; end
 xip = 2*xii/(1-eta)-1; etap = 2*eta-1;
%---

  Jc=0;
  for k=0:m
   for l=0:m-k
    Jc=Jc+1;
    pror(Jc) = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;
   end
  end

  psil = vani*pror';

  zp(i) = psil(node);

end
plot3 (xp,yp,zp,'k');
end

%---
% plot the triangle
%---

xx(1)=0.0; yy(1)=0.0; zz(1)=0.0;
xx(2)=1.0; yy(2)=0.0; zz(2)=0.0;
xx(3)=0.0; yy(3)=1.0; zz(3)=0.0;
xx(4)=0.0; yy(4)=0.0; zz(4)=0.0;

plot3(xx,yy,zz,'k');

%-----
% done
%-----

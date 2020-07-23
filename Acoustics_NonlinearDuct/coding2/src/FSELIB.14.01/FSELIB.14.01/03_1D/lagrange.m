clear all
close all

%==========
% Lagrange interpolation
% based on given data
%
% calls lagrange_interp
%==========

%-----
% choose
%-----
 
ibase = 3; % Cheb
ibase = 1; % evenly spaced
ibase = 2; % Lob

if(ibase==1)
 mv(1)=2;mv(2)=4;mv(3)=8;mv(4)=11;
 mv(1)=2;
 mv(2)=3;
 mv(3)=4;
 mv(4)=5;
 mv(5)=6;
 mv(6)=7;
 ms = 6;
end

if(ibase==2)
 mv(1)=2;
 mv(2)=3;
 mv(3)=4;
 mv(4)=5;
 mv(5)=6;
 mv(6)=7;
 ms = 6;
 axis([-1 1 -0.2 1])
end

%-----------------
% prepare
%-----------------

figure(1)
hold on
xlabel('\xi','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)
box on

nplot = 64;
for i=1:nplot+1
  xi(i) = -1.0+2.0*(i-1)/nplot;
  y(i) = 1.0/(1.0+25.0*xi(i)^2);
end
plot(xi,y,'r-','linewidth',3);

%============
for ipass=1:ms
%============

m = mv(ipass);

%---------------
% evenly spaced:
%---------------

if(ibase==1)

 for i=1:m+1
  xi(i) = -1.0+2.0*(i-1)/m;
 end

end

%--------
% lobatto:
%--------

if(ibase==2)

xi(1) = -1.0;
if(m>1)
  [zL, wL]= lobatto(m-1);
  for i=2:m
   xi(i) = zL(i-1);
  end
end
xi(m+1)=1.0;

end

%-----------
% Chebyshev:
%-----------

if(ibase==3)

 for i=1:m+1
  xi(i) = cos((i-1)*pi/m);
 end

end

%-----------------
% function values:
%-----------------

for i=1:m+1
  y(i) = 1.0/(1.0+xi(i)^2);
  y(i) = 1.0/(1.0+25.0*xi(i)^2);
end

a = xi(1);
b = xi(m+1);
nplot=2*64;
step = (b-a)/nplot;
 
for l=1:nplot+1
  xint(l) = a+step*(l-1);
  yint(l) = lagrange_interp(m,xi,y,xint(l));
end
 
if(ipass==1)
   plot(xint,yint,'k-','linewidth',1);
elseif(ipass==2)
   plot(xint,yint,'k--','linewidth',1);
elseif(ipass==3)
   plot(xint,yint,'k-.','linewidth',1);
elseif(ipass==4)
   plot(xint,yint,'k:','linewidth',1);
elseif(ipass==5)
   plot(xint,yint,'k--','linewidth',2);
elseif(ipass==6)
   plot(xint,yint,'k-.','linewidth',2);
elseif(ipass==7)
   plot(xint,yint,'r:','linewidth',5);
end

%=====
end % of ipass
%=====

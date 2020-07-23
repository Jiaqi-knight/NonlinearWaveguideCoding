%================================
% psi_modes
%
% prepare graphs of the element
% modal interpolation functions
% for a 3-node triangle
%===============================

close all
clear all
mode =  11;

mode = 112;
mode = 113;
mode = 123;
mode = 212;
mode = 213;
mode = 223;

%---
% prepare
%---

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
% plot the triangle
%---

figure(1)
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\zeta','fontsize',14);
set(gca,'fontsize',15)
box on
view([68,22])

xx(1)=0.0; yy(1)=0.0; zz(1)=0.0;
xx(2)=1.0; yy(2)=0.0; zz(2)=0.0;
xx(3)=0.0; yy(3)=1.0; zz(3)=0.0;
xx(4)=0.0; yy(4)=0.0; zz(4)=0.0;

plot3(xx,yy,zz,'k');

%---
% vertical lines
%---

for i=1:N
 for j=1:M+2-i

  xp(j) = x(i);
  yp(j) = y(j);  
  zeta= 1.0-xp(j)-yp(j);
  X = 2.0*xp(j)/(1.0001-yp(j))-1.0;      % xi'
  Y = 2.0*yp(j)-1.0;                     % eta'

      if(mode==112)
   zp(j)= xp(j)*zeta;
  elseif(mode==113)
   zp(j)= yp(j)*zeta;
  elseif(mode==123)
   zp(j)= xp(j)*yp(j);
  elseif(mode==11)
   zp(j)= xp(j)*yp(j)*zeta;      % bubble mode, m=3
  elseif(mode==212)
   zp(j)= 0.5*(1-X) * 0.5*(1+X)* (1-Y)^3/8.0 * 3*X;
  elseif(mode==213)
   zp(j)= 0.5*(1+Y) * 0.5*(1-X)*  0.5*(1-Y)  * 3*Y;
  elseif(mode==223)
   zp(j)= 0.5*(1+X) * 0.5*(1+Y)*  0.5*(1-Y)  * 3*Y;
  end
 
 end
 plot3(xp,yp,zp,'k');
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
  X = 2.0*xp(i)/(1-yp(i))-1.0; % xi'
  Y = 2.0*yp(i)-1.0;            % eta'

      if(mode==112)
   zp(i)= xp(i)*zeta; 
  elseif(mode==113)
   zp(i)= yp(i)*zeta;
  elseif(mode==123)
   zp(i)= xp(i)*yp(i);
  elseif(mode==11)
   zp(i)= xp(i)*yp(i)*zeta;    % bubble mode, m=3
  elseif(mode==212)
   zp(i)= 0.5*(1-X) *0.5*(1+X) *(1-Y)^3/8.0 * 3*X;
  elseif(mode==213)
   zp(i)= 0.5*(1+Y) *0.5*(1-X) *  0.5*(1-Y) * 3*Y;
  elseif(mode==223)
   zp(i)= 0.5*(1+X) *0.5*(1+Y) *  0.5*(1-Y) * 3*Y;
  end

 end
 plot3(xp,yp,zp,'k');
end
                                                                                
%-----
% done
%-----


clear all
close all

%=======
% modal interpolation functions
% on a quadrilateral
%=======

%---
% prepare
%---

figure(1)
hold on
xlabel('\xi','fontsize',15);
ylabel('\eta','fontsize',15);
zlabel('\zeta','fontsize',15);
set(gca,'fontsize',15)
view([36,18])

modex=2;
modey=2;

N=32;
M=32;

Dx=2.0/N;
Dy=2.0/M;

%-----
% lobatto nodes
%-----

xi(1)=-1; xi(2)=0; xi(3)=1;
et(1)=-1; et(2)=0; et(3)=1;

%-----
% mesh
%-----

for i=1:N+1 
   x(i)=-1.0+Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end
                                                                                
for j=1:M+1
   y(j)=-1.0+Dy*(j-1.0);
%   fprintf (1,'%6.2f',Y(j))
end

%-----
% graph
%-----

for i=1:N+1
for j=1:M+1
 if(modex==1 & modey==1)
  z(i,j) = (1-x(i))*(1+x(i))*(1-y(j))*(1+y(j))/16;
 elseif(modex==1 & modey==2)
  z(i,j) = (1-x(i))*(1+x(i))*(1-y(j))*(1+y(j))*3*y(j)/16;
 elseif(modex==2 & modey==1)
  z(i,j) = (1-x(i))*(1+x(i))*(1-y(j))*(1+y(j))*3*x(i)/16;
 elseif(modex==2 & modey==2)
  z(i,j) = (1-x(i))*(1+x(i))*(1-y(j))*(1+y(j))*9*x(i)*y(j)/16;
 end
end
end

mesh(x,y,z')   % plotting, labelling, and formatting

%---
% done
%---

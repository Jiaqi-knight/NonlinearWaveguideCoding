
clear all
close all

%=======
% psi_q22
%
% interpolation functions for the
% 2X2 quadrilateral
%=======

nodex = 2;
nodey = 1;

nodex = 1;
nodey = 1;

nodex = 1;
nodey = 2;

nodex = 2;
nodey = 2;

%---
% prepare
%---

figure(1)
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
view([40,18])
box on

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
     if(nodex==1 & nodey==1)
 z(i,j) = (x(i)-xi(2))*(x(i)-xi(3))/( (xi(1)-xi(2))*(xi(1)-xi(3)) ) ...
         *(y(j)-et(2))*(y(j)-et(3))/( (et(1)-et(2))*(et(1)-et(3)) );
 elseif(nodex==1 & nodey==2)
 z(i,j) = (x(i)-xi(2))*(x(i)-xi(3))/( (xi(1)-xi(2))*(xi(1)-xi(3)) ) ...
         *(y(j)-et(1))*(y(j)-et(3))/( (et(2)-et(1))*(et(2)-et(3)) );
 elseif(nodex==2 & nodey==1)
 z(i,j) = (x(i)-xi(1))*(x(i)-xi(3))/( (xi(2)-xi(1))*(xi(2)-xi(3)) ) ...
         *(y(j)-et(2))*(y(j)-et(3))/( (et(1)-et(2))*(et(1)-et(3)) );
 elseif(nodex==2 & nodey==2)
 z(i,j) = (x(i)-xi(1))*(x(i)-xi(3))/( (xi(2)-xi(1))*(xi(2)-xi(3)) ) ...
         *(y(j)-et(1))*(y(j)-et(3))/( (et(2)-et(1))*(et(2)-et(3)) );
 end
end
end

mesh(x,y,z')   % plotting, labelling, and formatting

%---
% done
%---

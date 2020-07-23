
close all
clear all

%=======
% psi_q12
%
% interpolation functions for the
% 12-node quadrilateral
%=======

figure(1)
hold on
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
view([-19 20])
box on

node = 1;
node = 2;

%---
% prepare
%---

N=32;
M=32;

Dx=2.0/N;
Dy=2.0/M;

for i=1:N+1 
   x(i)=-1.0+Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end
                                                                                
for j=1:M+1
   y(j)=-1.0+Dy*(j-1.0);
%   fprintf (1,'%6.2f',Y(j))
end

for i=1:N+1
for j=1:M+1
 if(node==1)
  z(i,j)= (1.0-x(i))*(1.0-y(j))*(-10+9.0*x(i)^2+9.0*y(j)^2)/32; 
 elseif(node==2)
  z(i,j)= 9.0*(1.0-x(i)^2)*(1.0-y(j))*(1.0-3.0*x(i))/32;
 end
end
end

mesh(x,y,z')   % plotting, labelling, and formatting


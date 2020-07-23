
clear all
close all

%=======
% psi_q8
%
% interpolation functions
% for the 8-node quadrilateral element
%=======

figure(1)
hold on
axis equal
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
zlabel('\psi','fontsize',14);
set(gca,'fontsize',14)
view([71,18])
box on

node = 2;
node = 1;


N=32;
M=32;

Dx=2.0/N;
Dy=2.0/M;

%------
% plotting mesh
%------

for i=1:N+1 
   x(i)=-1.0+Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end
                                                                                
for j=1:M+1
   y(j)=-1.0+Dy*(j-1.0);
%   fprintf (1,'%6.2f',Y(j))
end

%------
% graph
%------

for i=1:N+1
 for j=1:M+1
  if(node==1)
   z(i,j)= -0.25*(1.0-x(i))*(1.0-y(j))*(x(i)+y(j)+1.0);   % 8-node-corner
   axis([-1 1 -1 1 -0.5 1])
  elseif(node==2)
   z(i,j)= 0.5*(1.0-x(i)^2)*(1.0-y(j));   % 8-node mid side
   axis([-1 1 -1 1 -0.0 1.1])
  end
 end
end

mesh(x,y,z)   % plotting, labelling, and formatting

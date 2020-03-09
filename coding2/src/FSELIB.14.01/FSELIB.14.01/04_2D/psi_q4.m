
clear all
close all

%=======================================
% Code phi_q4
%
% element-node interpolation functions
% over the canonical rectangular element
%=======================================

%-----------------
% prepare the grid
%-----------------

N = 32;
M = 32;

Dx = 2.0/N;
Dy = 2.0/M;

for i=1:N+1 
   x(i) =-1.0+Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end
                                                                                
for j=1:M+1
   y(j) =-1.0+Dy*(j-1.0);
%   fprintf (1,'%6.2f',Y(j))
end

%-----
% plot
%-----

for i=1:N+1
 for j=1:M+1
  z(i,j) = 0.25*(1.0-x(i))*(1-y(j));   % 4-node
 end
end

%---
% plot
%---

figure(1)
mesh(x,y,z)   % plotting, labelling, and formatting
box on
axis equal
set(gca,'fontsize',15)
axis equal
xlabel('\xi','fontsize',15);
ylabel('\eta','fontsize',15);
zlabel('\psi_1','fontsize',15);

%-----
% done
%-----

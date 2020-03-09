
clear all
close all

%================
% Code lapl6_d_rc
%
% Solution of Laplace's equation
% in a rectangle with a circular hole
% with the Dirichlet boundary condition
%
% a: hole radius
%================

%-----
% input data
%-----

ndiv = 2;
a = 0.5;
NQ = 6;    % gauss-triangle quadrature

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl6_rc(a, ndiv);

%-------------------------------------
% specify Dirichlet boundary conditions
%-------------------------------------

for i=1:ng
 if(gfl(i)==1)     % outer boundary
   bcd(i) = 1.0;   % example
 end
 if(gfl(i)==2)     % inner boundary
   bcd(i) = 0.0;   % example
 end
end

%-------------------------------------
% assemble the global diffusion matrix
%-------------------------------------

gdm = zeros(ng,ng); % initialize

area = 0.0;

for l=1:ne          % loop over the elements

% compute the element diffusion matrix

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[edm_elm, arel] = edm6 ...
...
     (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
     ,NQ);

area = area + arel;

   for i=1:6
     i1 = c(l,i);
     for j=1:6
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
   end

end

% disp (area)

%------------------------
% set the right-hand side
% and implement the Dirichlet BC
%------------------------

for i=1:ng
 b(i) = 0.0;
end

for j=1:ng
 if(gfl(j) ~= 0) 
   for i=1:ng 
    b(i) = b(i) - gdm(i,j) * bcd(j);
    gdm(i,j) = 0; gdm(j,i) = 0;
   end
%   disp (j) 
   gdm(j,j) = 1.0; b(j) = bcd(j);
 end

end

% disp (ng)
% disp (gdm)
% disp (b)

%------------------------
% solve the linear system
%------------------------

f = b/gdm'; 

%  disp (f)

%-----
% plot
%-----

plot_6 (ne,ng,p,c,f);
axis([-1.1 3.1 -2.1 2.1]);

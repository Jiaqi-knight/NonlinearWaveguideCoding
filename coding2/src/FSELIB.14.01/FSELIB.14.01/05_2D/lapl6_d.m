
clear all
close all

%=============
% Code lapl6_d
%
% FSELIB
%
% Finite element code for solving
% Laplace's equation
% in a disk-like domain subject to
% the Dirichlet boundary condition,
% using 6-node triangles.
%=============

%-----------
% prepare
%-----------

figure(1)
hold on;
axis square
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('f','fontsize',14);
set(gca,'fontsize',15)
%axis([-1 1 -1 1]);
%axis([-0.8 0.8 -0.8 0.6]);
view([-7 18])
box on

%-----------
% input data
%-----------

ndiv = 2;  % triangulation level

NQ = 6;    % gauss-triangle quadrature

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl6_disk(ndiv);

%--------
% deform
%--------

defx = 0.60;

for i=1:ng
  p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
end

%-----------------------------------------------
% specify the Dirichlet boundary condition (dbc)
%-----------------------------------------------

for i=1:ng
 dbc(i) = 0.0; % default
 if(gfl(i)==1)
   dbc(i) = 1.0;            % another example
   dbc(i) = sin(pi*p(i,2));    % example
   dbc(i) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
   dbc(i) = p(i,1);            % linear in x
 end
end

%-------------------------------------
% assemble the global diffusion matrix
% and compute the domain surface area (optional)
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

area = area+arel;

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
 if(gfl(j)==1) 
   for i=1:ng 
    b(i) = b(i) - gdm(i,j) * dbc(j);
    gdm(i,j) = 0; gdm(j,i) = 0;
   end
%   disp (j) 
   gdm(j,j) = 1.0; b(j) = dbc(j);
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
  [f' p(:,1) gfl' dbc']

%-----
% plot
%-----

plot_6 (ne,ng,p,c,f);

%-----
% generate the extended connectivity matrix
% for 3-node sub-triangles
%-----

Ic=0;

for i=1:ne
 Ic=Ic+1;
 c3(Ic,1)=c(i,1); c3(Ic,2)=c(i,4); c3(Ic,3)=c(i,6);
 Ic=Ic+1;
 c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,2); c3(Ic,3)=c(i,5);
 Ic=Ic+1;
 c3(Ic,1)=c(i,5); c3(Ic,2)=c(i,3); c3(Ic,3)=c(i,6);
 Ic=Ic+1;
 c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,5); c3(Ic,3)=c(i,6);
end

trimesh (c3,p(:,1),p(:,2),f);

%-----
% done
%-----

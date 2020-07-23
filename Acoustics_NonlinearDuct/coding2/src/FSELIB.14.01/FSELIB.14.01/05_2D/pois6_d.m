
close all
clear all

%=============
% Code pois6_d
%
% Finite element code for solving Poisson's equation
% in a disk-like domain subject to
% the Dirichlet boundary condition,
% using 6-node triangles.
%=============

%-----------
% input data
%-----------

ndiv = 3;  % triangulation level

NQ = 6;    % gauss-triangle quadrature

ishape = 1;
ishape = 2;
ishape = 4;
ishape = 5;
ishape = 6;
ishape = 3;

%------------
% triangulate
%------------

if(ishape==1)
 ndiv = 3; 
 [ne,ng,p,c,efl,gfl] = trgl6_disk(ndiv);
elseif(ishape==2)
 ndiv = 3;
 [ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);
elseif(ishape==3)
 ndiv = 3;
 [ne,ng,p,c,efl,gfl] = trgl6_L(ndiv);
elseif(ishape==4)
 ndiv = 2;
 a = 0.5;
 [ne,ng,p,c,efl,gfl] = trgl6_sc(a, ndiv);
elseif(ishape==5)
 ndiv = 2;
 a = 0.5;
 [ne,ng,p,c,efl,gfl] = trgl6_ss(a,ndiv);
elseif(ishape==6)
 ndiv = 2;
 a = 0.5;
 [ne,ng,p,c,efl,gfl] = trgl6_rc (a,ndiv)
end

%--------
% deform
%--------

defx = 0.60; defy = 0.10;
defx = 0.00; defy = 0.00;

for i=1:ng
  p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
  p(i,2) = p(i,2)*(1.0-defy*p(i,2)^2 );
end

%----------------
% deform: azimuthal waves
%----------------
                                                                                
amp = 0.20;
amp = 0.00;
                                                                                
nphi=3;
                                                                                
for i=1:ng
 rad  = sqrt(p(i,1)^2+p(i,2)^2);
 if(rad<0.00001)
  angle = 0.10;
 else
  angle = acos(p(i,1)/rad);
 end
 radd = 1.0+amp*cos(nphi*angle);
 p(i,1) = p(i,1)*radd;
 p(i,2) = p(i,2)*radd;
end

%-----------------------------------------
% specify the Dirichlet boundary condition
%-----------------------------------------

for i=1:ng
 if(gfl(i)==1)
   bcd(i) = p(i,1);            % another example
   bcd(i) = sin(pi*p(i,2));    % example
   bcd(i) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
   bcd(i) = 1.0;            % another example
   bcd(i) = 0.0;            % another example
 end
 if(gfl(i)==2)
   bcd(i,2) = 0.0;            % another example
 end
end

%-------------------
% specify the source
%-------------------
                                                                                
for i=1:ng
  source(i) =-1.0;
end

%-------------------------------------
% assemble the global diffusion matrix
% and global mass matrix
% and compute the domain surface area (optional)
%-------------------------------------

gdm = zeros(ng,ng); % initialize
gmm = zeros(ng,ng); % initialize

area = 0.0;

for l=1:ne          % loop over the elements

% compute the element diffusion and mass matrices

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[edm_elm, emm_elm, arel] = edmm6 ...
...
     (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
     ,NQ);

area = area+arel;

   for i=1:6
     i1 = c(l,i);
     for j=1:6
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
   end

end

% disp (area)

%------------------------
% set the right-hand side
%------------------------

for i=1:ng
 b(i) = 0.0;
 for j=1:ng
  b(i) = b(i)-gmm(i,j)*source(j);
 end
end

%------------------------
% implement the Dirichlet BC
%------------------------

for j=1:ng
 if(gfl(j)==1|gfl(j)==2) 
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

%------------------
% plot the solution
%------------------

figure(1)
hold on
view([-19,18])
%axis equal
box on
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('f','fontsize',14);
set(gca,'fontsize',14)

%axis([-0.8 0.8 -0.8 0.6]);

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
%trisurf (c3,p(:,1),p(:,2),f,f);

%-----------------------------
% average the element gradient
% to compute the gradient
% at the global nodes
%-----------------------------

for j=1:ng              % initialize
  ggradx(j)=0.0;
  ggrady(j)=0.0;
 gitally(j)=0;
end

%---
% loop over the element nodes and compute
% the nodal gradient
%---
                                                                                
for l=1:ne
                                                                                
j=c(l,1); x(1)=p(j,1); y(1)=p(j,2); fn(1)=f(j);
j=c(l,2); x(2)=p(j,1); y(2)=p(j,2); fn(2)=f(j);
j=c(l,3); x(3)=p(j,1); y(3)=p(j,2); fn(3)=f(j);
j=c(l,4); x(4)=p(j,1); y(4)=p(j,2); fn(4)=f(j);
j=c(l,5); x(5)=p(j,1); y(5)=p(j,2); fn(5)=f(j);
j=c(l,6); x(6)=p(j,1); y(6)=p(j,2); fn(6)=f(j);

[gradx, grady] = elm6_grad (x,y,fn);

 for k=1:6
   j=c(l,k);
   gitally(j) = gitally(j)+1;
   ggradx(j) = ggradx(j)+gradx(k);
   ggrady(j) = ggrady(j)+grady(k);
 end

end

for j=1:ng
 ggradx(j)=ggradx(j)/gitally(j);
 ggrady(j)=ggrady(j)/gitally(j);
 ggradn(j) = sqrt(ggradx(j)^2+ggrady(j)^2);
end

figure(2)
plot_6 (ne,ng,p,c,ggradx);
figure(3)
plot_6 (ne,ng,p,c,ggrady);
figure(4)
plot_6 (ne,ng,p,c,ggradn);

%-----
% done
%-----

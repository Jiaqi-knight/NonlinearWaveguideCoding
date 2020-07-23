%===========
% Code psa6
%
% Plane stress analysis in a square
% or a square with a hole using 6-node triangles
%
% dependencies: elm6_abc elm_6_interp esm6 gauss_trgl psa6_stress
%=============

clear all
close all

%-----------
% prepare
%-----------

figure(1)
hold on
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
set(gca,'fontsize',14)
axis equal
box on
%axis([-1.2 1.2 -1.2 1.2]);
%axis([-24.2 24.2 -24.2 24.2]);

%-----------
% input data
%-----------

E = 1.0;  % modulus of elasticity
%E = 30000;  % modulus of elasticity
nu = 0.25; % Poisson ratio
nu = 0.50; % Poisson ratio
NQ = 7;   % gauss-triangle quadrature
ndiv = 2; % discretization level

ishape = 1;
ishape = 2;

%------------
% triangulate
%------------

if(ishape==1)

 [ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);   % square

elseif(ishape==2)

 a= 0.5;  % hole radius
 [ne,ng,p,c,efl,gfl] = trgl6_sc(a, ndiv); % square with a hole

end

disp('Number of elements:'); ne

%--------
% prepare
%--------

ng2 = 2*ng;
nus = nu^2;  % square of the Poisson ratio

%-------------------------
% count the interior nodes
%-------------------------

Ic = 0;

for j=1:ng
 if(gfl(j,1)==0) Ic=Ic+1; end
end

disp('Number of interior nodes:'); Ic

%--------------------------------------------
% specify the displacement boundary condition
% on the top, bottom, and left side, and 
% the traction on the right side
%--------------------------------------------

for i=1:ng

  gfl(i,2) = 0;  % initialize

 if(gfl(i,1)==1)

  % default:  zero traction BC

  gfl(i,2) = 1;  gfl(i,3) = 0.0; % fx = 0
                 gfl(i,4) = 0.0; % fy = 0

%-----
% apply a load at right edge
%-----

  if(p(i,1) > 0.99 )
    gfl(i,3) = 0.0;                 % fx
    gfl(i,4) = 0.01*(1.0-p(i,2)^2); % fy
%   gfl(i,4) = 40*(1.0-p(i,2)^2);   % fy
  end

%-----
% apply a load at top
%-----

%  if(p(i,2) > 0.99 )
%    gfl(i,3) = 0.0;                  % fx
%    gfl(i,4) = -0.10*(1.0-p(i,1)^2); % fy
%  end

%-----
% zero displacement at the left edge
%-----

  if(p(i,1) < -0.99)  
    gfl(i,2) = 2;    % displacement BC
    gfl(i,3) = 0.0;  % (vx) x-displacement 
    gfl(i,4) = 0.0;  % (vy) y-displacement
  end

%-----
% zero displacement at right edge
%-----

%  if(p(i,1) > 0.99)  
%    gfl(i,2) = 2;    % displacement BC
%    gfl(i,3) = 0.0;  % (vx) x-displacement 
%    gfl(i,4) = 0.0;  % (vy) y-displacement
%  end

%-----
% zero displacement at bottom
%-----

%  if(p(i,2) < -0.99)  
%    gfl(i,2) = 2;    % displacement BC
%    gfl(i,3) = 0.0;  % (vx) x-displacement 
%    gfl(i,4) = 0.0;  % (vy) y-displacement
%  end

 end

end

%----------------------
% deform to a rectangle
%----------------------

defx =24.0; defy=6.0;  % deformation factors
defx = 2.0; defy=1.0;  % deformation factors
defx = 0.25;defy=1.0;  % deformation factors
defx = 4.0; defy=1.0;  % deformation factors
defx = 2.0; defy=1.0;  % deformation factors
defx = 1.0; defy=1.0;  % deformation factors

for i=1:ng
  p(i,1) = p(i,1)*defx;
  p(i,2) = p(i,2)*defy;
end

%------------------
% plot the elements
%------------------

for i=1:ne

 i1=c(i,1); i2=c(i,2); i3=c(i,3); i4=c(i,4); i5=c(i,5); i6=c(i,6);

 xp(1)=p(i1,1); xp(2)=p(i4,1); xp(3)=p(i2,1); xp(4)=p(i5,1);
                xp(5)=p(i3,1); xp(6)=p(i6,1); xp(7)=p(i1,1);
 yp(1)=p(i1,2); yp(2)=p(i4,2); yp(3)=p(i2,2); yp(4)=p(i5,2);
                yp(5)=p(i3,2); yp(6)=p(i6,2); yp(7)=p(i1,2);
% plot(xp, yp);
 plot(xp, yp,'k+','markersize',5)
end

%---------------------------------------
% assemble the global stiffness martices
%---------------------------------------

gsm_xx = zeros(ng,ng); % initialize
gsm_xy = zeros(ng,ng); % initialize
gsm_yy = zeros(ng,ng); % initialize

for l=1:ne          % loop over the elements

% compute the element stiffness matrices

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[esm_xx, esm_xy, esm_yy, arel] = esm6 ...
...
     (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ,NQ);

   for i=1:6
     i1 = c(l,i);
     for j=1:6
       j1 = c(l,j);
       gsm_xx(i1,j1) = gsm_xx(i1,j1) + esm_xx(i,j);
       gsm_xy(i1,j1) = gsm_xy(i1,j1) + esm_xy(i,j);
       gsm_yy(i1,j1) = gsm_yy(i1,j1) + esm_yy(i,j);
     end
   end

end

%-------------------------
% form the transpose block
%-------------------------

gsm_yx = gsm_xy';    % a prime denotes the transpose

%--------------------------
% define the block matrices
%--------------------------

Rxu =    gsm_xx/(1-nus) + 0.5*gsm_yy/(1+nu);
Rxv = nu*gsm_xy/(1-nus) + 0.5*gsm_yx/(1+nu);
Ryu = Rxv';
Ryv =    gsm_yy/(1-nus) + 0.5*gsm_xx/(1+nu);

%-------------------------
% compile the grand matrix
%-------------------------

for i=1:ng
 for j=1:ng
  Gm(i,   j) = E*Rxu(i,j); Gm(i,   ng+j) = E*Rxv(i,j);
  Gm(ng+i,j) = E*Ryu(i,j); Gm(ng+i,ng+j) = E*Ryv(i,j);
 end
end

% disp (Gm);

%-----------------------------
% assemble the right-hand side
% of the linear system
%-----------------------------

for i=1:ng2
 b(i) = 0.0;
end

for i=1:ne

   % first side

  l1=c(i,1); l2=c(i,2); l3=c(i,4);

  if(gfl(l1,2)==1 & gfl(l2,2)==1 )
   side = sqrt( (p(l2,1)-p(l1,1))^2+(p(l2,2)-p(l1,2))^2 );
   b(l1) = b(l1) + gfl(l1,3)*side*1.0/6.0;
   b(l2) = b(l2) + gfl(l2,3)*side*1.0/6.0;
   b(l3) = b(l3) + gfl(l3,3)*side*2.0/3.0;
   k1 = l1+ng;
   k2 = l2+ng;
   k3 = l3+ng; 
   b(k1) = b(k1) + gfl(k1,4)*side*1.0/6.0;
   b(k2) = b(k2) + gfl(k2,4)*side*1.0/6.0;
   b(k3) = b(k3) + gfl(k3,4)*side*2.0/3.0;
  end

   % second  side

  l1=c(i,2); l2=c(i,3); l3=c(i,5);

  if(gfl(l1,2)==1 & gfl(l2,2)==1 )
   side = sqrt( (p(l2,1)-p(l1,1))^2+(p(l2,2)-p(l1,2))^2 );
   b(l1) = b(l1) + gfl(l1,3)*side*1.0/6.0;
   b(l2) = b(l2) + gfl(l2,3)*side*1.0/6.0;
   b(l3) = b(l3) + gfl(l3,3)*side*2.0/3.0;
   k1 = l1+ng; k2 = l2+ng; k3 = l3+ng; 
   b(k1) = b(k1) + gfl(l1,4)*side*1.0/6.0;
   b(k2) = b(k2) + gfl(l2,4)*side*1.0/6.0;
   b(k3) = b(k3) + gfl(l3,4)*side*2.0/3.0;
  end

   % third side

  l1=c(i,3); l2=c(i,1); l3=c(i,6);

  if(gfl(l1,2)==1 & gfl(l2,2)==1 )
   side = sqrt( (p(l2,1)-p(l1,1))^2+(p(l2,2)-p(l1,2))^2 );
   b(l1) = b(l1) + gfl(l1,3)*side*1.0/6.0;
   b(l2) = b(l2) + gfl(l2,3)*side*1.0/6.0;
   b(l3) = b(l3) + gfl(l3,3)*side*2.0/3.0;
   k1 = l1+ng; k2 = l2+ng; k3 = l3+ng; 
   b(k1) = b(k1) + gfl(k1,4)*side*1.0/6.0;
   b(k2) = b(k2) + gfl(k2,4)*side*1.0/6.0;
   b(k3) = b(k3) + gfl(k3,4)*side*2.0/3.0;
  end

end

%---------------------------
% implement the Dirichlet BC
%---------------------------

for j=1:ng  % run over global nodes

%---
% if(gfl(j,2)==1)           % boundary node with traction BC
%  b(j)    = gfl(j,3);
%  b(j+ng) = gfl(j,4);
% end
%---

%---
 if(gfl(j,2)==2)  % boundary node with displacement BC

   for i=1:ng2
    b(i) = b(i) - Gm(i,j) * gfl(j,3) - Gm(i,ng+j) * gfl(j,4);
    Gm(i,j) = 0; Gm(i,ng+j) = 0;
    Gm(j,i) = 0; Gm(ng+j,i) = 0;
   end

  Gm(j,  j)     = 1.0;
  Gm(j+ng,j+ng) = 1.0;

  b(j)    = gfl(j,3);
  b(j+ng) = gfl(j,4);

 end
%---

end

%disp (Gm);

%------------------------
% solve the linear system
%------------------------

f = b/Gm'; 

%-------------------------
% assign the displacements
%-------------------------

for i=1:ng
 vx(i) = f(i);
 vy(i) = f(ng+i);
end

%---------------------------
% plot the deformed elements
%---------------------------

for i=1:ng
  p(i,1) = p(i,1)+vx(i);
  p(i,2) = p(i,2)+vy(i);
end

for i=1:ne

 i1=c(i,1); i2=c(i,2); i3=c(i,3);i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xp(1)=p(i1,1); xp(2)=p(i4,1); xp(3)=p(i2,1); xp(4)=p(i5,1);
                xp(5)=p(i3,1); xp(6)=p(i6,1); xp(7)=p(i1,1);
 yp(1)=p(i1,2); yp(2)=p(i4,2); yp(3)=p(i2,2); yp(4)=p(i5,2);
                yp(5)=p(i3,2); yp(6)=p(i6,2); yp(7)=p(i1,2);
 plot(xp, yp,'-ko','markersize',5)
end

%-----
% reset
%-----

for i=1:ng
  p(i,1) = p(i,1)-vx(i);
  p(i,2) = p(i,2)-vy(i);
end

%-----------------------------
% average the element nodal stresses
% to compute the stresses
% at the global nodes
%-----------------------------

for j=1:ng              % initialize the global node stresses
 gsig_xx(j) = 0.0;
 gsig_xy(j) = 0.0;
 gsig_yy(j) = 0.0;
 gitally(j) = 0;
end

%---
% loop over the element nodes and compute
% the element nodal stresses
%---

for l=1:ne

j=c(l,1); xe(1)=p(j,1); ye(1)=p(j,2); vxe(1)=vx(j); vye(1)=vy(j);
j=c(l,2); xe(2)=p(j,1); ye(2)=p(j,2); vxe(2)=vx(j); vye(2)=vy(j);
j=c(l,3); xe(3)=p(j,1); ye(3)=p(j,2); vxe(3)=vx(j); vye(3)=vy(j);
j=c(l,4); xe(4)=p(j,1); ye(4)=p(j,2); vxe(4)=vx(j); vye(4)=vy(j);
j=c(l,5); xe(5)=p(j,1); ye(5)=p(j,2); vxe(5)=vx(j); vye(5)=vy(j);
j=c(l,6); xe(6)=p(j,1); ye(6)=p(j,2); vxe(6)=vx(j); vye(6)=vy(j);

[sig_xx, sig_xy, sig_yy] = psa6_stress ...
...
     (xe,ye,vxe,vye,E,nu);

 for k=1:6
   j=c(l,k);
   gitally(j) = gitally(j) + 1;
   gsig_xx(j) = gsig_xx(j) + sig_xx(k);
   gsig_xy(j) = gsig_xy(j) + sig_xy(k);
   gsig_yy(j) = gsig_yy(j) + sig_yy(k);
 end

end

for j=1:ng
 gsig_xx(j) = gsig_xx(j)/gitally(j);
 gsig_xy(j) = gsig_xy(j)/gitally(j);
 gsig_yy(j) = gsig_yy(j)/gitally(j);
end

%-----
% generalize the connectivity matrix
% for 3-node sub-triangles
%-----
                                                                                
Ic=0;

for i=1:ne
 Ic=Ic+1;
 c3(Ic,1) = c(i,1); c3(Ic,2) = c(i,4); c3(Ic,3) = c(i,6);
 Ic=Ic+1;
 c3(Ic,1) = c(i,4); c3(Ic,2) = c(i,2); c3(Ic,3) = c(i,5);
 Ic=Ic+1;
 c3(Ic,1) = c(i,5); c3(Ic,2) = c(i,3); c3(Ic,3) = c(i,6);
 Ic=Ic+1;
 c3(Ic,1) = c(i,4); c3(Ic,2) = c(i,5); c3(Ic,3) = c(i,6);
end

%----------------------
% plot the stress field
% using the matlab trimesh function
%----------------------

figure(2)
hold on
trimesh (c3,p(:,1),p(:,2),gsig_xy);
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('\sigma_{xy}','fontsize',14);
set(gca,'fontsize',14)
box on
axis([-1.1 1.1 -1.1 1.1 -0.005 0.015])
axis([-1.1 1.1 -1.1 1.1 -0.010 0.025])
view([33 44])

figure(3)
hold on
trimesh (c3,p(:,1),p(:,2),gsig_xx);
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('\sigma_{xx}','fontsize',14);
set(gca,'fontsize',14)
box on
axis([-1.1 1.1 -1.1 1.1 -0.04 0.060])
axis([-1.1 1.1 -1.1 1.1 -0.04 0.080])
view([33 44])

figure(4)
hold on
trimesh (c3,p(:,1),p(:,2),gsig_yy);
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('\sigma_{yy}','fontsize',14);
set(gca,'fontsize',14)
box on
axis([-1.1 1.1 -1.1 1.1 -0.015 0.020])
axis([-1.1 1.1 -1.1 1.1 -0.025 0.040])
view([33 44])

%-----------------------------------
% plot the displacement vector field
%-----------------------------------

figure(77)
hold on
axis equal
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
set(gca,'fontsize',14)
%axis([-1.2 1.2 -1.2 1.2]);
%axis([-24.2 24.2 -24.2 24.2]);

%-----
% done
%-----

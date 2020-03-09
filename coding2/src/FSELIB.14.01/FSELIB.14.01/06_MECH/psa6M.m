%===========
% CODE psa6M
%
% In-plane deformation of a membrane patch under
% a uniform body force in plane stress analysis
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
nu = 0.5; % Poisson ratio
NQ = 6;   % gauss-triangle quadrature
bx =  1.00;
by =  0.00;

ishape = 2;
ishape = 1;
ishape = 3;

%------------
% triangulate
%------------

if(ishape==1)
 ndiv = 2; % discretization level
 [ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);   % square
elseif(ishape==2)
 a = 0.5;  % hole radius
 ndiv = 2; % discretization level
 [ne,ng,p,c,efl,gfl] = trgl6_sc(a, ndiv); % square with a hole
elseif(ishape==3)
  ndiv = 3; % discretization level
 [ne,ng,p,c,efl,gfl] = trgl6_disk(ndiv)
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

ic = 0;

for j=1:ng
 if(gfl(j,1)==0) ic=ic+1; end
end

disp('Number of interior nodes:'); ic

%-------------------------------------------------
% specify the zero displacement boundary condition
%-------------------------------------------------

for i=1:ng
 if(gfl(i,1)==1)    % boundary node
   gfl(i,2) = 2;    % displacement BC
   gfl(i,3) = 0.0;  % u x-displacement 
   gfl(i,4) = 0.0;  % v y-displacement
 end
end

%-------
% deform
%-------

def = 0.5;   % deformation factor
def = 0.0;   % deformation factor

for i=1:ng
  p(i,1)=p(i,1)*(1.0 + def);
  p(i,2)=p(i,2)*(1.0 - def);
end

%------------------
% plot the elements
%------------------

for i=1:ne

 i1=c(i,1);i2=c(i,2); i3=c(i,3);i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xp(1)=p(i1,1); xp(2)=p(i4,1); xp(3)=p(i2,1); xp(4)=p(i5,1);
                xp(5)=p(i3,1); xp(6)=p(i6,1); xp(7)=p(i1,1);
 yp(1)=p(i1,2); yp(2)=p(i4,2); yp(3)=p(i2,2); yp(4)=p(i5,2);
                yp(5)=p(i3,2); yp(6)=p(i6,2); yp(7)=p(i1,2);
  plot(xp, yp,'k+','markersize',5)
end

%------------------------------
% assemble the global stiffness
% and matrices
%------------------------------

gsm_xx = zeros(ng,ng); % initialize
gsm_xy = zeros(ng,ng); % initialize
gsm_yy = zeros(ng,ng); % initialize

for i=1:ng
  gmass(i) = 0.0;
end

for l=1:ne          % loop over the elements

% compute the element stiffness matrices
% and RHS

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[esm_xx, esm_xy, esm_yy, emass, arel] = esmm6 ...
...
     (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ);

   for i=1:6
     i1 = c(l,i);
     gmass(i1) = gmass(i1) + emass(i);
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
  Gm(i,   j) = E*Rxu(i,j); Gm(i,   j+ng) = E*Rxv(i,j);
  Gm(ng+i,j) = E*Ryu(i,j); Gm(ng+i,j+ng) = E*Ryv(i,j);
 end
end

% disp (Gm);

%------------------------
% set the right-hand side of the linear system
% and implement the Dirichlet BC
%------------------------

for i=1:ng2
 b(i) = 0.0;
end

for j=1:ng  % run over nodes

b(j)    = gmass(j)*bx;
b(j+ng) = gmass(j)*by;

%---
 if(gfl(j,1)== 1 & gfl(j,2)==2)  % boundary node with displacement BC

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

%---
% reset
%---

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
 gsig_xx(j)=0.0;
 gsig_xy(j)=0.0;
 gsig_yy(j)=0.0;
 gitally(j)=0;
end

%---
% loop over the element nodes and compute
% the element nodal stresses
%---

for l=1:ne

j=c(l,1); xe(1)=p(j,1); ye(1)=p(j,2); ue(1)=vx(j); ve(1)=vy(j);
j=c(l,2); xe(2)=p(j,1); ye(2)=p(j,2); ue(2)=vx(j); ve(2)=vy(j);
j=c(l,3); xe(3)=p(j,1); ye(3)=p(j,2); ue(3)=vx(j); ve(3)=vy(j);
j=c(l,4); xe(4)=p(j,1); ye(4)=p(j,2); ue(4)=vx(j); ve(4)=vy(j);
j=c(l,5); xe(5)=p(j,1); ye(5)=p(j,2); ue(5)=vx(j); ve(5)=vy(j);
j=c(l,6); xe(6)=p(j,1); ye(6)=p(j,2); ue(6)=vx(j); ve(6)=vy(j);

[sig_xx, sig_xy, sig_yy] = psa6_stress ...
...
     (xe,ye,ue,ve,E,nu);

 for k=1:6
   j=c(l,k);
   gitally(j) = gitally(j) + 1;
   gsig_xx(j) = gsig_xx(j) + sig_xx(k);
   gsig_xy(j) = gsig_xy(j) + sig_xy(k);
   gsig_yy(j) = gsig_yy(j) + sig_yy(k);
 end

end

for j=1:ng
 gsig_xx(j)=gsig_xx(j)/gitally(j);
 gsig_xy(j)=gsig_xy(j)/gitally(j);
 gsig_yy(j)=gsig_yy(j)/gitally(j);
end

%-----
% extend the connectivity matrix
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
trimesh (c3,p(:,1),p(:,2),gsig_xx);
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('\sigma_{xx}','fontsize',14);
set(gca,'fontsize',14)
box on
axis([-1.1 1.1 -1.1 1.1 -0.80 0.80])
view([33 44])

figure(3)
hold on
trimesh (c3,p(:,1),p(:,2),gsig_xy);
quiver (p(:,1)',p(:,2)',vx,vy);
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('\sigma_{xy}','fontsize',14)
set(gca,'fontsize',14)
box on
axis([-1.1 1.1 -1.1 1.1 -0.40 0.40])
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
axis([-1.1 1.1 -1.1 1.1 -0.40 0.40])
view([33 44])


%-----
% done
%-----

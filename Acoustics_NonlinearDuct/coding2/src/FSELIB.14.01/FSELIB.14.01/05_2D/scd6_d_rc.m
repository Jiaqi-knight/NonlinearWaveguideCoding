
close all
clear all

%=============================================
% Code scd6_d_rc
%
% FSELIB
%
% Steady convection--diffusion equation
% with the Dirichlet boundary condition
% in a rectangular domain with a circular hole
% using 6-node triangles
%=============================================

%-----
% input data
%-----

k =1.0; rho=1.0; cp=1.0;
NQ = 6;
ndiv = 2;
a = 0.25;  % cylinder radius

U = 20.0;  % Pe = aU/kappa
U = 00.0;  % Pe = aU/kappa
U = 10.0;  % Pe = aU/kappa

%----------
% constants
%----------

kappa = k/(rho*cp);

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl6_rc(a, ndiv);

%-----------------------------------------
% specify Dirichlet the boundary condition
%-----------------------------------------

for i=1:ng

 if(gfl(i)==1)  % outer boundary
   bcd(i) = 1.0;
 end

 if(gfl(i)==2)  % inner boundary
   bcd(i) = 0.0;
 end

end

%---------------------------------------------------
% assemble the global diffusion and advection matrix
%---------------------------------------------------

gdm = zeros(ng,ng); % initialize
gam = zeros(ng,ng); % initialize

for l=1:ne % Loop over the elements

% compute the element diffusion 
% and advection matrices

  j=c(l,1); x1=p(j,1); y1=p(j,2);
  j=c(l,2); x2=p(j,1); y2=p(j,2);
  j=c(l,3); x3=p(j,1); y3=p(j,2);
  j=c(l,4); x4=p(j,1); y4=p(j,2);
  j=c(l,5); x5=p(j,1); y5=p(j,2);
  j=c(l,6); x6=p(j,1); y6=p(j,2);

% Define the node velocities:

  [u1, v1] = scd6_vel(U,a,x1,y1);
  [u2, v2] = scd6_vel(U,a,x2,y2);
  [u3, v3] = scd6_vel(U,a,x3,y3);
  [u4, v4] = scd6_vel(U,a,x4,y4);
  [u5, v5] = scd6_vel(U,a,x5,y5);
  [u6, v6] = scd6_vel(U,a,x6,y6);

[edm_elm, eam_elm, arel] = edam6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6 ...
   ,u1,v1, u2,v2, u3,v3, u4,v4, u5,v5, u6,v6 ...
   ,NQ);

   for i=1:6
     i1 = c(l,i);
     for j=1:6
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gam(i1,j1) = gam(i1,j1) + eam_elm(i,j);
     end
   end
end

%-------------------
% coefficient matrix
%-------------------

lsm = kappa*gdm+gam;

%------------------------
% set the right-hand side and
% implement the boundary conditions
%------------------------

for i=1:ng
 b(i) = 0.0;
end

for j=1:ng

 if(gfl(j) == 1 | gfl(j) == 2) 

   for i=1:ng 
    b(i) = b(i) - lsm(i,j) * bcd(j);
    lsm(i,j) = 0; lsm(j,i) = 0;
   end

   lsm(j,j) = 1.0;
   b(j) = bcd(j);

 end

end

%------------------------
% solve the linear system
%------------------------

f = b/lsm'; 

%---------
% plotting
%---------

figure(1)
hold on
view([-19,18])
axis equal
box on
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('f','fontsize',14);
set(gca,'fontsize',14)
axis([-1 3 -1 1 0 1.2])
view([29 20])

plot_6 (ne,ng,p,c,f);

%-----------------------------
% extended connectivity matrix
% for 3-node sub-triangles
%-----------------------------

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
                                                                                
%---------
% trimesh plotting
%---------

trimesh (c3,p(:,1),p(:,2),f);

ido = 0;

if(ido==1)

 for j=1:ng
  xpp(1) = p(j,1); ypp(1)=p(j,2); 
  if(gfl(j) == 0) plot(xpp, ypp,'o');end
  if(gfl(j) == 1) plot(xpp, ypp,'x');end
  if(gfl(j) == 2) plot(xpp, ypp,'+');end
 end

end

%----
% done
%---

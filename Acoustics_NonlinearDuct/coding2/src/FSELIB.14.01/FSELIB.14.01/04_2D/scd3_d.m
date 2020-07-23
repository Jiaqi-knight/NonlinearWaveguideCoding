
clear all
close all

%======================================================
% Code scd3_d
%
% Finite-element code for steady convection--diffusion
% subject to the Dirichlet boundary condition
%======================================================

%-----------
% input data
%-----------

k = 1.0; rho=1.0; cp=1.0;
U = 00.0;
U = 02.0;
U = 05.0;
U = 10.0;
V = 00.0;
ndiv = 3;

%----------
% constants
%----------

kappa = k/(rho*cp);

%-----------------------
% triangulate and deform
%-----------------------

 [ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);
%[ne,ng,p,c,efl,gfl] = trgl3_delaunay;

%-------
% deform
%-------

defx = 0.5;
defx = 0.0;

for i=1:ng
 p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
end

%-----------------------------------------
% specify the Dirichlet boundary condition
%-----------------------------------------

for i=1:ng
 if(gfl(i)==1)
%  bcd(i) = 1.0;  % example
   bcd(i) = sin(pi*p(i,2));    % example
   bcd(i) = p(i,1)*sin(0.5*pi*p(i,2));    % example
   bcd(i) = p(i,1);  % example
 end
end

%------------------------------
% assemble the global diffusion
% and advection matrix
%------------------------------

gdm = zeros(ng,ng); % initialize
gam = zeros(ng,ng); % initialize

for l=1:ne % Loop over the elements

% compute the element diffusion 
% and advection matrices

  j=c(l,1); x1=p(j,1); y1=p(j,2);
  j=c(l,2); x2=p(j,1); y2=p(j,2);
  j=c(l,3); x3=p(j,1); y3=p(j,2);

  [edm_elm] = edm3 (x1,y1,x2,y2,x3,y3);

  [eam_elm] = eam3 (U,V,x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
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

%----------------------------
% compute the right-hand side
%----------------------------

for i=1:ng
 b(i) = 0.0;
end

for m=1:ng
 if(gfl(m)==1) 
   for i=1:ng 
    b(i) = b(i) - lsm(i,m) * bcd(m);
    lsm(i,m) = 0; lsm(m,i) = 0;
   end
   lsm(m,m) = 1.0; b(m) = bcd(m);
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
plot_3 (ne,ng,p,c,f);
trimesh (c,p(:,1),p(:,2),f);
%trisurf (c,p(:,1),p(:,2),f,f);
view([-19,18])
box on
%axis([-0.8 0.8 -0.8 0.6]);

%-----
% done
%-----

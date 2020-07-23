close all
clear all

%========================================
% Code lapl3_eig
%
% Eigenfunctions of Laplace's equation
% in a disk-like domain
% using 3-node triangles
%========================================

%-----------
% input data
%-----------

ndiv = 1;  % discretization level
ndiv = 2;  % discretization level
ndiv = 3;  % discretization level
ndiv = 5;  % discretization level
ndiv = 4;  % discretization level

%------------
% triangulate
%------------

 [ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);
%[ne,ng,p,c,efl,gfl] = trgl3_delaunay;

%-------
% deform
%-------

defx = 0.6;
defx = 0.0;

for i=1:ng
 p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
end

%-------------------------------------
% assemble the global diffusion and mass matrices
%-------------------------------------

gdm = zeros(ng,ng); % initialize
gmm = zeros(ng,ng); % initialize

for l=1:ne          % loop over the elements

% compute the element diffusion and mass matrices

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);

[edm_elm] = edm3 (x1,y1,x2,y2,x3,y3);
[emm_elm] = emm3 (x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
   end

end

%--------------------
% reduce the matrices by removing equations
% corresponding to boundary nodes
%--------------------

Ic=0;

for i=1:ng

 if(gfl(i)==0)
  Ic=Ic+1;
  map(Ic) = i;

  Jc=0;
  for j=1:ng
   if(gfl(j)==0)
    Jc=Jc+1;
    A(Ic,Jc) = gdm(i,j);
    B(Ic,Jc) = gmm(i,j);
   end
  end

 end

end

ngred = Ic;

%------------------------
% compute the eigenvalues
%------------------------

[V,D] = eig(A,B);

[2.405^2 5.520^2]

%-------------
% map and plot
%-------------

for order=1:6

D(order,order)

 for i=1:ng
  f(i) = 0.0;
 end

 for i=1:ngred
 f(map(i)) = V(i,order);
 end

 figure(order)
 plot_3 (ne,ng,p,c,f);
 trimesh (c,p(:,1),p(:,2),f);
 %trisurf (c,p(:,1),p(:,2),f,f);
 view([-19,18])
 box on
 %axis([-0.8 0.8 -0.8 0.6]);

end

%-----
% done
%-----


close all
clear all

%========================================
% Code lapl6_eig
%
% Eigenfunctions of Laplace's equation
% using 6-node triangles
%========================================

%-----------
% input data
%-----------

ndiv = 5;  % discretization level
ndiv = 4;  % discretization level
ndiv = 0;  % discretization level
ndiv = 1;  % discretization level
ndiv = 3;  % discretization level
ndiv = 2;  % discretization level

NQ = 12;    % gauss-triangle quadrature

ishape = 1;
ishape = 4;

%------------
% triangulate
%------------

if(ishape==1)
 [ne,ng,p,c,efl,gfl] = trgl6_disk (ndiv);
elseif(ishape==2)
 a=0.5;
 [ne,ng,p,c,efl,gfl] = trgl6_ss(a, ndiv);
elseif(ishape==3)
 a=0.5;
 [ne,ng,p,c,efl,gfl] = trgl6_sc(a, ndiv);
elseif(ishape==4)
 [ne,ng,p,c,efl,gfl] = trgl6_L(ndiv)
end

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
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[edm_elm, emm_elm, arel] = edmm6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ);

   for i=1:6
     i1 = c(l,i);
     for j=1:6
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

%-----
% map and plot
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

for order=1:6

D(order,order)

 for i=1:ng
  f(i) = 0.0;
 end

 for i=1:ngred
  f(map(i)) = V(i,order);
 end

 figure(order)
 plot_6 (ne,ng,p,c,f);
 trimesh (c3,p(:,1),p(:,2),f);
 view([-19,18])
 view([60,24])
 box on
 %axis([-0.8 0.8 -0.8 0.6]);

end

%-----
% done
%-----

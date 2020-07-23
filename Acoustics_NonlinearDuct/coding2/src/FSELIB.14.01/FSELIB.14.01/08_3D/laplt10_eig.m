close all
clear all

%========================================
% Code laplt10_eig
%
% Eigenfunctions of Laplace's equation
% in a spherical-like domain
% using 10-node tetrahedra
%========================================

%-----------
% input data
%-----------

ndiv = 3;  % discretization level
ndiv = 0;  % discretization level
ndiv = 2;  % discretization level
ndiv = 1;  % discretization level

%-----------
% discretize
%-----------

[ne,ng,p,c,efl,gfl] = tetra10_sphere8 (ndiv);

%-------
% deform
%-------

defx = 0.00;
defy = 0.00;
defz = 0.00;

for i=1:ng
 p(i,1) = p(i,1)*(1.0+defx*p(i,1)^2 );
 p(i,2) = p(i,2)*(1.0+defy*p(i,2)^2 );
 p(i,3) = p(i,3)*(1.0+defz*p(i,3)^2 );
end

%-------------------------------------
% assemble the global diffusion and mass matrices
%-------------------------------------

gdm = zeros(ng,ng); % initialize
gmm = zeros(ng,ng); % initialize

volume = 0;

for l=1:ne          % loop over the elements

% compute the element diffusion and mass matrices

j=c(l,1); x1 =p(j,1); y1 =p(j,2); z1 =p(j,3);
j=c(l,2); x2 =p(j,1); y2 =p(j,2); z2 =p(j,3);
j=c(l,3); x3 =p(j,1); y3 =p(j,2); z3 =p(j,3);
j=c(l,4); x4 =p(j,1); y4 =p(j,2); z4 =p(j,3);
j=c(l,5); x5 =p(j,1); y5 =p(j,2); z5 =p(j,3);
j=c(l,6); x6 =p(j,1); y6 =p(j,2); z6 =p(j,3);
j=c(l,7); x7 =p(j,1); y7 =p(j,2); z7 =p(j,3);
j=c(l,8); x8 =p(j,1); y8 =p(j,2); z8 =p(j,3);
j=c(l,9); x9 =p(j,1); y9 =p(j,2); z9 =p(j,3);
j=c(l,10);x10=p(j,1); y10=p(j,2); z10=p(j,3);

[edm_elm, emm_elm, volel] = edmm_t10 ...
...
   (x1,y1,z1, x2,y2,z2, x3,y3,z3 ...
   ,x4,y4,z4, x5,y5,z5, x6,y6,z6 ...
   ,x7,y7,z7, x8,y8,z8, x9,y9,z9, x10,y10,z10);

   for i=1:10
     i1 = c(l,i);
     for j=1:10
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
   end

   volume = volume + volel;

end

%--------------------
% clip the matrices by removing equations
% corresponding to boundary nodes
%--------------------

Ic=0;

for i=1:ng

 if(gfl(i)==0)

  Ic = Ic+1;
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

eigenvalues = eig(A,B);
eigenvalues(1:5)

[pi^2] % smallest eigenvalue

%-----
% done
%-----

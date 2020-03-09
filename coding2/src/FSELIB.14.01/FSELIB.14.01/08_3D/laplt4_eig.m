%========================================
% Code laplt4_eig
%
% Eigenfunctions of Laplace's equation
% in a spherical-like domain
% using 4-node tetrahedra
%========================================

close all
clear all

%-----------
% input data
%-----------

ishape = 1;
ishape = 2;

%-----------
% discretize
%-----------

if(ishape==1)
 ndiv = 2;  % discretization level
 ndiv = 3;  % discretization level
 ndiv = 0;  % discretization level
 ndiv = 1;  % discretization level
 [ne,ng,p,c,efl,gfl] = tetra4_sphere8 (ndiv);
elseif(ishape==2)
 ndiv = 0;  % discretization level
 ndiv = 1;  % discretization level
 ndiv = 2;  % discretization level
 ndiv = 3;  % discretization level
 [ne,ng,p,c,efl,gfl] = tetra4_cube8 (ndiv);
end

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

for l=1:ne          % loop over the elements

% compute the element diffusion and mass matrices

j=c(l,1); x1=p(j,1); y1=p(j,2); z1=p(j,3);
j=c(l,2); x2=p(j,1); y2=p(j,2); z2=p(j,3);
j=c(l,3); x3=p(j,1); y3=p(j,2); z3=p(j,3);
j=c(l,4); x4=p(j,1); y4=p(j,2); z4=p(j,3);

[edm_elm, vlm_elm] = edm_t4 (x1,y1,z1,x2,y2,z2 ...
                            ,x3,y3,z3,x4,y4,z4);

[emm_elm] = emm_t4  (x1,y1,z1,x2,y2,z2 ...
                     ,x3,y3,z3,x4,y4,z4);

   for i=1:4
     i1 = c(l,i);
     for j=1:4
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
   end

end

%--------------------
% clip the matrices by removing equations
% corresponding to boundary nodes
%--------------------

Ic=0;

for i=1:ng

 if(gfl(i)==0)

  Ic=Ic+1;
  map(Ic) = i;

  Jc = 0;
  for j=1:ng
   if(gfl(j)==0)
    Jc = Jc+1;
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

if(ishape==1)
 [pi^2] % smallest eigenvalue
end

%-----
% done
%-----

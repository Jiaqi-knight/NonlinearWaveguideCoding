%================================================
% Code lapl3_dn_sqr
%
% Finite element solution of Laplace's equation
% in a deformed square
% with Dirichlet and Neumann boundary conditions
%===============================================

close all
clear all

%-----------
% input data
%-----------

k = 1.0;      % conductivity;
ndiv = 3;     % level of triangulation

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl3_sqr(ndiv);

%---------------------------------------
% find the element sides on the boundary 
%---------------------------------------

Ic = 0;   % counter

for i=1:ne   % run over the elements

 efl(i,4) = efl(i,1);  % wrap around the element
 c(i,4) = c(i,1);

 for j=1:3
  if( efl(i,j)==1 & efl(i,j+1)==1 )
   Ic = Ic+1;
   face1(Ic) = c(i,j);
   face2(Ic) = c(i,j+1);
  end
 end

end

nbe = Ic;    % number of boundary sides

%---------------------------------------
% specify the Neumann boundary condition
% on the left side 
%---------------------------------------

for i=1:nbe  % run along the sides 

 m = face1(i);
 l = face2(i);
 face3(i) = 1;

 if(abs(p(m,1)-1)<0.00001 | abs(p(l,1)-1)< 0.00001 )   % Neumann side
   face3(i)= 2;     % Neumann flag
   bcn1(i) = 2.0;   % Neumann condition at first  side point
   bcn2(i) = 2.0;   % Neumann condition at second side point
 end

end

%-----------------------------------------
% specify the Dirichlet boundary condition
% on the upper, right, and lower sides
%-----------------------------------------

for i=1:ng      % initialize node flag
 bcf(i) = 0;    % bcf(i)=1 will indicate a Dirichlet node
end

for i=1:nbe

 m = face1(i);
 l = face2(i);

 if(abs(p(m,2)-1) <0.00001 | abs(p(l,2)-1)<0.00001 )   % Dirichlet side
   bcf(m)=1;
   bcf(l)=1;    % Dirichlet flag
   bcd(m)=0.0;  % Dirichlet condition at first  side point
   bcd(l)=0.0;  % Dirichlet condition at second side point
 end

 if(abs(p(m,2)+1) <0.00001 | abs(p(l,2)+1)<0.00001 )   % Dirichlet side
   bcf(m)=1;
   bcf(l)=1;    % Dirichlet flag
   bcd(m)=0.0;  % Dirichlet condition at first  side point
   bcd(l)=0.0;  % Dirichlet condition at second side point
 end

 if(abs(p(m,1)+1)<0.00001 | abs(p(l,1)+1)< 0.00001 )   % Dirichlet side
   bcf(m)=1;
   bcf(l)=1;    % Dirichlet flag
   bcd(m)=0.0;  % Dirichlet condition at first  side point
   bcd(l)=0.0;  % Dirichlet condition at second side point
 end

end

%-----------------------------------------
% deform the square into a frivolous shape
%-----------------------------------------

for i=1:ng
 rad = sqrt(p(i,1)^2+p(i,2)^2);
 p(i,1)=p(i,1)*(0.25+0.00*rad);
 p(i,2)=p(i,2)*(1.00-0.10*rad);
end

%-------------------------------------
% assemble the global diffusion matrix
%-------------------------------------

gdm = zeros(ng,ng); % initialize

for l=1:ne    % loop over the elements

% compute the element matrices

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);

[edm_elm] = edm3(x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
   end
end

%-------------------
% initialize the rhs
%-------------------

for i=1:ng
 b(i) = 0.0;
end

%----------------------------------------
% Neumann integral on the right-hand side
%-----------------------------------------

for i=1:nbe

 if(face3(i)==2)

  m = face1(i);
  l = face2(i);
  xe1 = p(m,1); ye1 = p(m,2);
  xe2 = p(l,1); ye2 = p(l,2);
  edge = sqrt((xe2-xe1)^2+(ye2-ye1)^2);
  int1 = edge * ( bcn1(i)/3 + bcn2(i)/6 );
  int2 = edge * ( bcn1(i)/6 + bcn2(i)/3 );
  b(m) = b(m)+int1/k;
  b(l) = b(l)+int2/k;

 end

end

%-------------------------------------------
% implement the Dirichlet boundary condition
%-------------------------------------------

for j=1:ng

 if(bcf(j)==1)

   for i=1:ng
    b(i) = b(i) - gdm(i,j) * bcd(j);
    gdm(i,j) = 0.0; gdm(j,i) = 0.0;
   end
   gdm(j,j) = 1.0; b(j) = bcd(j);

 end

end

%------------------------
% solve the linear system
%------------------------

f = b/gdm'; 

%---------
% plotting
%---------

figure(1)
plot_3 (ne,ng,p,c,f);
trimesh (c,p(:,1),p(:,2),f);
%trisurf (c,p(:,1),p(:,2),f,f);
view([-19,18])
box on
axis equal
axis([-0.3 0.3 -1.1 1.1 0 0.8]);

%-----
% done
%-----

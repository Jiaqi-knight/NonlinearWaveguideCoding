%=======================================
% Code lapl3_dn
%
% Solution of Laplace's equation in a
% disk-like domain with Dirichlet
% and Neumann boundary conditions
%======================================

close all
clear all

%-----------
% input data
%-----------

k = 1.0;      % conductivity
ndiv = 3;     % level of triangulation

%-----------------------
% triangulate and deform
%-----------------------

[ne,ng,p,c,efl,gfl] = trgl3_disk(ndiv);

%-------
% deform
%-------

defx = 0.00;
defx = 0.30;

for i=1:ng
 p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
end

%-----------------------------------------
% find the element edges on the boundary C
%-----------------------------------------

Ic = 0;   % counter

for i=1:ne

 efl(i,4) = efl(i,1);     % wrap around the element
   c(i,4) =   c(i,1);     % wrap around the element

 for j=1:3
  if( efl(i,j)==1 & efl(i,j+1)==1 )
   Ic = Ic+1;
   face1(Ic) = c(i,j);
   face2(Ic) = c(i,j+1);
  end
 end

end

nbe = Ic    % number of boundary edges

%---------------------------------------
% specify the Neumann boundary condition
% on the right side (example)
%---------------------------------------

for i=1:nbe  % run along the boundary edges

 m = face1(i);
 l = face2(i);
     face3(i)=1; % default Dirichlet flag

 if(p(m,1)>0.00001 | p(l,1)> 0.00001 )   % Neumann edge
%  if(p(m,2)>0.00001 | p(l,2)> 0.00001 )   % Neumann edge (quadrant)
   face3(i)= 2;    % Neumann flag
   bcn1(i)= 1.0;   % Neumann condition at first  edge node
   bcn2(i)= 1.0;   % Neumann condition at second edge node
%  end
 end

end

%-----------------------------------------
% specify the Dirichlet boundary condition
% on the left half of the solution domain
%-----------------------------------------

for i=1:ng      % initialize node flag
 bcf(i) = 0;    % bcf(i)=1 will indicate a Dirichlet node
end

for i=1:nbe

 m = face1(i);
 l = face2(i);
 if(p(m,1)<-0.00001 | p(l,1)<-0.00001 )   % Dirichlet side
   bcf(m)=1;   % Dirichlet flag
   bcf(l)=1;   % Dirichlet flag
   bcd(m)=0.0;           % Dirichlet condition at first  edge node
   bcd(l)=0.0;           % Dirichlet condition at second edge node
 end

% quadrant
%
% if(p(m,2)<-0.00001 | p(l,2)<-0.00001 )   % Dirichlet side
%   bcf(m)=1; bcf(l)=1; % Dirichlet flag
%   bcd(m)=0.0;           % Dirichlet condition at first  edge node
%   bcd(l)=0.0;           % Dirichlet condition at second edge node
% end

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
    gdm(i,j) = 0; gdm(j,i) = 0;
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
axis([-1.1 1.1 -1.1 1.1 0 2.0]);

%-----
% done
%-----

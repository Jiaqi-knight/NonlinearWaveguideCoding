%====================================
% Code hlm3_n
%
% Code for Helmholtz's equation in a
% disk-like domain with the Neumann
% boundary condition
%====================================

close all
clear all

%-----------
% input data
%-----------

k = 1.0;      % conductivity;
ndiv = 3;     % level of triangulation

alpha =16.0;  % Helmholtz coefficient;
alpha =05.0;  % Helmholtz coefficient;
alpha =01.0;  % Helmholtz coefficient;

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl3_disk(ndiv);

%-------
% deform
%-------

defx = 0.30;
defx = 0.00;

for i=1:ng
 p(i,1)=p(i,1)*(1.0-defx*p(i,2)^2 );
end

%---------------------------------------
% specify the Neumann boundary condition
%---------------------------------------

for i=1:ng
 if(gfl(i)==1)
  bcn(i) = 1.0;  % constant
 end
end

%----------------------------------------------
% assemble the global diffusion and mass matrix
%----------------------------------------------

gdm = zeros(ng,ng); % initialize

gmm = zeros(ng,ng); % initialize

for l=1:ne    % loop over the elements

% compute the element matrices

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);

[edm_elm] = edm3(x1,y1,x2,y2,x3,y3);
[emm_elm] = emm3(x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
       gmm(i1,j1) = gmm(i1,j1) + emm_elm(i,j);
     end
   end

end

%---------------------------
% set up the right-hand side
%---------------------------

for i=1:ng
 b(i) = 0.0;  % initialize
end

%----------------------------------------
% Neumann integral on the right-hand side
%----------------------------------------

for i=1:ne

 efl(i,4) = efl(i,1); c(i,4) = c(i,1);  % replicate the first node

 for j=1:3         % run around the sides

  if( efl(i,j)==1 & efl(i,j+1)==1 )

    j1=c(i,j);
    j2=c(i,j+1);

    xe1 = p(j1,1); ye1 = p(j1,2);
    xe2 = p(j2,1); ye2 = p(j2,2);

    edge = sqrt((xe2-xe1)^2+(ye2-ye1)^2);

    int1 = edge * ( bcn(j1)/3 + bcn(j2)/6 );
    int2 = edge * ( bcn(j1)/6 + bcn(j2)/3 );

    b(j1) = b(j1)+int1/k;
    b(j2) = b(j2)+int2/k;

  end

 end
end

%-------------------
% coefficient matrix
%-------------------

lsm = gdm-alpha*gmm; 

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

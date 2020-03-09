%========================================
% Code pois3_d
%
% Solution of Poisson's equation
% with the Dirichlet boundary condition
% in a disk-like domain
% using 3-node triangles
%========================================

%-----------
% input data
%-----------

ndiv = 4;  % discretization level

%------------
% triangulate
%------------

%[ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);
[ne,ng,p,c,efl,gfl] = trgl3_sqr(ndiv);

%-------
% deform
%-------

defx = 0.60;
defx = 0.00;
defy = 0.00;

for i=1:ng
 p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
 p(i,2) = p(i,2)*(1.0-defy*p(i,2)^2 );
end

%----------------
% azimuthal waves
%----------------

amp = 0.20;

nphi=3;

for i=1:ng
 rad  = sqrt(p(i,1)^2+p(i,2)^2);
 if(rad<0.00001)
  angle = 0.10;
 else
  angle = acos(p(i,1)/rad);
 end
 radd = 1.0+amp*cos(nphi*angle);
 p(i,1) = p(i,1)*radd;
 p(i,2) = p(i,2)*radd;
end

%-----------------------------------------
% specify the Dirichlet boundary condition
%-----------------------------------------

for i=1:ng
 if(gfl(i,1)==1)
   gfl(i,2) = sin(pi*p(i,2));    % example
   gfl(i,2) = p(i,1);            % another example
   gfl(i,2) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
   gfl(i,2) = 0.0;   % another example
 end
end

%-------------------
% specify the source
%-------------------

for i=1:ng
  source(i) =-10.0;
end

%-------------------------------------
% assemble the global diffusion matrix
% and global mass matrix
%-------------------------------------

gdm = zeros(ng,ng); % initialize
gmm = zeros(ng,ng); % initialize

for l=1:ne          % loop over the elements

% compute the element diffusion 
% and mass matrix

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

% disp (gdm);
% disp (gmm);

%---------------------------------------------
% set the right-hand side of the linear system
% and implement the Dirichlet boundry condition
%----------------------------------------------

for i=1:ng
 b(i) = 0.0;
 for j=1:ng
  b(i) = b(i)-gmm(i,j)*source(j);
 end
end

for j=1:ng
 if(gfl(j,1)==1) 
   for i=1:ng 
    b(i) = b(i) - gdm(i,j)*gfl(j,2);
    gdm(i,j)=0; gdm(j,i)=0;
   end
   gdm(j,j)=1.0;
   b(j)=gfl(j,2);
 end
end

%------------------------
% solve the linear system
%------------------------

f = b/gdm'; 

%-----
% plot
%-----

plot_3 (ne,ng,p,c,f);
hold on
trimesh (c,p(:,1),p(:,2),f);
%trisurf (c,p(:,1),p(:,2),f,f);

%-----
% done
%-----

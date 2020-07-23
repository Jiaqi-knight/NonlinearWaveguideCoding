close all
clear all

%========================================
% Code lapl3_d
%
% Solution of Laplace's equation
% with the Dirichlet boundary condition
% in a disk-like domain
% using 3-node triangles
%========================================

%-----------
% prepare
%-----------

figure(1)
hold on;
axis square
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
zlabel('f','fontsize',14);
set(gca,'fontsize',15)
%axis([-1 1 -1 1]);
%axis([-0.8 0.8 -0.8 0.6]);
view([-19,18])
box on

%-----------
% input data
%-----------

ndiv = 2;  % discretization level
ndiv = 0;  % discretization level
ndiv = 1;  % discretization level

%------------
% triangulate
%------------

 [ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);
%[ne,ng,p,c,efl,gfl] = trgl3_delaunay;

%-------
% deform
%-------

defx = 0.0;
defx = 0.6;

for i=1:ng
 p(i,1) = p(i,1)*(1.0-defx*p(i,2)^2 );
 p(i,1) = p(i,1)*(1.0+defx);
end

%-----------------------------------------
% specify the Dirichlet boundary condition (bcd)
%-----------------------------------------

for i=1:ng
 if(gfl(i)==1)
   bcd(i) = sin(pi*p(i,2));    % example
   bcd(i) = p(i,1)^2;            % another example
   bcd(i) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
   bcd(i) = p(i,1);            % linear in x
 end
end

%-------------------------------------
% assemble the global diffusion matrix
%-------------------------------------

gdm = zeros(ng,ng); % initialize

for l=1:ne          % loop over the elements

% compute the element diffusion matrix

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);

[edm_elm] = edm3 (x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
   end

end

% disp (gdm);

%---------------------------------------------
% set the right-hand side of the linear system
% and implement the Dirichlet boundry condition
%----------------------------------------------

for i=1:ng
 b(i) = 0.0;
end

for j=1:ng

 if(gfl(j)==1) 
   for i=1:ng 
    b(i) = b(i) - gdm(i,j)*bcd(j);
    gdm(i,j) = 0; gdm(j,i) = 0;
   end
   gdm(j,j)=1.0;
   b(j) = bcd(j);
 end

end

%------------------------
% solve the linear system
%------------------------

f = b/gdm'; 

[f' p(:,1) gfl' bcd']

%-----
% plot
%-----

plot_3 (ne,ng,p,c,f);
trimesh (c,p(:,1),p(:,2),f);
%trisurf (c,p(:,1),p(:,2),f,f);

%-----
% done
%-----

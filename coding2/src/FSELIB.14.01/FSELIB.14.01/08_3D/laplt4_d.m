close all
clear all

%========================================
% Code lapl4_d
%
% Solution of Laplace's equation with
% the Dirichlet boundary condition
% using 4-node tetrahedral elements
%========================================

%--------
% prepare
%--------

 figure(1)
 hold on
 axis equal
 xlabel('x','fontsize',14);
 ylabel('y','fontsize',14);
 zlabel('z','fontsize',14);
 set(gca,'fontsize',14)
 view([-35 18])
 box on
 axis off
%axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

%-----------
% input data
%-----------
                                                                                
ndiv = 3;  % discretization level
ndiv = 0;  % discretization level
ndiv = 1;  % discretization level
ndiv = 2;  % discretization level

ishape = 2; % cube
ishape = 1; % sphere

%-----------
% discretize
%-----------

if(ishape==1)
% [ne,ng,p,c,efl,gfl] = tetra4_sphere12 (ndiv);
  [ne,ng,p,c,efl,gfl] = tetra4_sphere8 (ndiv);
elseif(ishape==2)
  [ne,ng,p,c,efl,gfl] = tetra4_cube8 (ndiv);
end

%-------
% deform
%-------

defx = 0.20;
defy = 0.20;
defz = 0.20;

for i=1:ng
 p(i,1) = p(i,1)*(1.0+defx*p(i,1)^2 );
 p(i,2) = p(i,2)*(1.0+defy*p(i,2)^2 );
 p(i,3) = p(i,3)*(1.0+defz*p(i,3)^2 );
end

%-----------------------------------------------
% specify the Dirichlet boundary condition (bcd)
%-----------------------------------------------

for i=1:ng
 bcd(i) = 0.0;   % default
 if(gfl(i)==1)
   bcd(i) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
   bcd(i) = sin(pi*p(i,2));    % example
   bcd(i) = 0.0;            % constant
   bcd(i) = 1.4;            % constant
   bcd(i) = p(i,2);            % linear in y
   bcd(i) = sin(pi*(p(i,1)+p(i,2)));    % example
   bcd(i) = p(i,1);            % linear in x
 end
end
                                                                                
%-------------------------------------
% assemble the global diffusion matrix
%-------------------------------------
                                                                                
gdm = zeros(ng,ng); % initialize

volume = 0.0;

for l=1:ne          % loop over the elements

% compute the element diffusion matrix
                                                                                
j=c(l,1); x1=p(j,1); y1=p(j,2); z1=p(j,3);
j=c(l,2); x2=p(j,1); y2=p(j,2); z2=p(j,3);
j=c(l,3); x3=p(j,1); y3=p(j,2); z3=p(j,3);
j=c(l,4); x4=p(j,1); y4=p(j,2); z4=p(j,3);

[edm_elm, vlm_elm] = edm_t4 (x1,y1,z1,x2,y2,z2 ...
                            ,x3,y3,z3,x4,y4,z4);

%edm_elm
%pause

  volume = volume+vlm_elm;

 %vlm_elm
 %pause

   for i=1:4
     i1 = c(l,i);
     for j=1:4
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
   end

end

volume
% disp(gdm);
%gdm

%-----------------------------------------------
% set the right-hand side of the linear system
% and implement the Dirichlet boundary condition
%-----------------------------------------------

for i=1:ng
 rhs(i) = 0.0;
end
                                                                                
for j=1:ng
 if(gfl(j)==1)
   for i=1:ng
     rhs(i) = rhs(i) - gdm(i,j)*bcd(j);
     gdm(i,j) = 0.0; gdm(j,i) = 0.0;
   end
   gdm(j,j) = 1.0;
   rhs(j) = bcd(j);
 end
end

%------------------------
% solve the linear system
%------------------------

f = rhs/gdm';

%-----
% plot
%-----

plot_t4 (ne,ng,p,c,f)

%---
% test
%---

itest = 0;
itest = 1;

if(itest==1)

 for i=1:ng
  ftest(i)= p(i,1);
 end

 res = gdm*ftest'-rhs';

 [res f' p(:,1) gfl' bcd']

 figure(55)
 hold on
 plot(p(:,1), f, 'ro')
 plot([-1 1],[-1 1],'k--')
 box on

end

%-----
% done
%-----

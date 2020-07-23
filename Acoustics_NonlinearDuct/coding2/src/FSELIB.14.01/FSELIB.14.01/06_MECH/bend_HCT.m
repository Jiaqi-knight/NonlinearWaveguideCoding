clear all
close all

%==========================
% Code bend_HCT
%
% Solution of the biharmonic equation
% for plate bending
% using the HCT element
%==========================

%-----
% input data
%-----

E_B = 1.0; % bending modulus
NQ = 4;      % quadrature order
w0 = 4;    % load

ndiv = 3;  % discretization level
ndiv = 2;  % discretization level

%---
% prepare
%---

figure(1)
hold on
%axis square
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('f','fontsize',14)
set(gca,'fontsize',14)
box on
view([23,26])

%------------
% triangulate
%------------

% [ne,ng,p,c,efl,gfl] = trgl6_disk(ndiv);
  [ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);

%-------
% deform
%-------

defx = 2.0;
defx = 1.0;

for i=1:ng
  p(i,1) = defx*p(i,1);
end

%--------------------------------------------
% make sure the mid-nodes are along the edges
%--------------------------------------------

for i=1:ne
 for j=1:2
  p(c(i,4),j) = 0.5*( p(c(i,1),j)+p(c(i,2),j) );
  p(c(i,5),j) = 0.5*( p(c(i,2),j)+p(c(i,3),j) );
  p(c(i,6),j) = 0.5*( p(c(i,3),j)+p(c(i,1),j) );
 end
end

%--------------------------------------
% mark the global vertex and edge nodes
%--------------------------------------

for i=1:ne
   gfln(c(i,1)) = 0; gfln(c(i,2)) = 0; gfln(c(i,3)) = 0;  % vertex nodes
   gfln(c(i,4)) = 1; gfln(c(i,5)) = 1; gfln(c(i,6)) = 1;  % edge nodes
end

%--------------------------------------
% count the number of global edge nodes
%--------------------------------------

nge=0;  % number of edge nodes

for i=1:ng
 if(gfln(i) == 1)
   nge = nge+1;
   xpp(nge)=p(i,1); ypp(nge)=p(i,2);
 end
end

plot(xpp, ypp,'^');

ic=0;
for i=1:ng
 if(gfln(i) == 0)
   ic=ic+1; xppp(ic)=p(i,1); yppp(ic)=p(i,2);
 end
end

plot(xppp, yppp,'ko');

%-----------------------------
% total number of global modes
%-----------------------------

ngv = ng-nge;     % vertex nodes

ngm = 3*ngv+nge;  % number of global modes

%------------------------------------------
% mark the position of the rows and columns
% of the modes corresponding to the global nodes
%------------------------------------------

for i=1:ng
  s(i) = 1;
  for j=1:i-1
   if(gfln(j) == 1)    % edge node, 1 dof
     s(i) = s(i)+1;
   else                % vertex node, 3 dofs
     s(i) = s(i)+3; 
   end
  end
end

%-----
% initialize the orientation index of the normal derivative
% at the midside nodes
%-----

for i=1:ng
 Idid(i)=1.0;
end

%---------------------------------------------
% assemble the global bending stiffness matrix
% and rhs
%---------------------------------------------

gbsm = zeros(ngm,ngm); % initialize
b    = zeros(1,ngm);    % right-hand side

area = 0;

for l=1:ne          % loop over the elements

% compute the element bending stiffness matrix and rhs

j1=c(l,1); j2=c(l,2); j3=c(l,3);   % global element labels
j4=c(l,4); j5=c(l,5); j6=c(l,6);

x1=p(j1,1); y1=p(j1,2);   % three vertices
x2=p(j2,1); y2=p(j2,2);
x3=p(j3,1); y3=p(j3,2);

[ebsm_elm, rhs, arel] = HCT_ebsm (x1,y1, x2,y2, x3,y3, NQ, w0);

e(1) = s(j1); e(2) = s(j1)+1; e(3) = s(j1)+2;   % position in the global matrix
e(4) = s(j2); e(5) = s(j2)+1; e(6) = s(j2)+2;
e(7) = s(j3); e(8) = s(j3)+1; e(9) = s(j3)+2;
e(10)= s(j4); e(11)= s(j5);   e(12)= s(j6);

   for i=1:12

     fci = 1.0;
         if(i==10) fci = Idid(j4);
     elseif(i==11) fci = Idid(j5);
     elseif(i==12) fci = Idid(j6); end

     for j=1:12
      fcj = 1.0;
          if(j==10) fcj = Idid(j4);
      elseif(j==11) fcj = Idid(j5);
      elseif(j==12) fcj = Idid(j6); end
       gbsm(e(i),e(j)) = gbsm(e(i),e(j)) + fci*fcj*ebsm_elm(i,j);
     end
     b(e(i)) = b(e(i)) - fci*rhs(i);
   end

 Idid(j4) = -1; Idid(j5) = -1; Idid(j6) = -1;

area = area+arel;

end

% disp (gbsm);

%---------------------------------------
% Implement the homogeneous Dirichlet BC
% No need to modify the RHS
%---------------------------------------

for j=1:ng

 if(gfl(j,1)==1)   % boundary node 

   j1 = s(j);
   j2 = s(j)+1;
   j3 = s(j)+2;

   for i=1:ngm

    gbsm(i,j1)=0.0; gbsm(j1,i)=0.0;

    if(gfln(j) == 0)   % vertex node
      gbsm(i,j2)=0.0; gbsm(j2,i)=0.0;
      gbsm(i,j3)=0.0; gbsm(j3,i)=0.0;
    end
   end

   gbsm(j1,j1) = 1.0; b(j1) = 0.0;

   if(gfln(j) == 0)   % vertex node
     gbsm(j2,j2) = 1.0; b(j2) = 0.0;
     gbsm(j3,j3) = 1.0; b(j3) = 0.0;
   end

 end
end

%------------------------
% solve the linear system
%------------------------

sol = b/gbsm'; 

%---------------------------------------
% extract the deflection at the vertices
%---------------------------------------

for i=1:ng
  if(gfln(i) == 0)
    f(i) = sol(s(i));
  end
end

%-------------------------------
% interpolate the midside values
%-------------------------------

for l=1:ne

 j1=c(l,1); j2=c(l,2); j3=c(l,3);   % global element labels
 j4=c(l,4); j5=c(l,5); j6=c(l,6);

 e(1) = s(j1); e(2) = s(j1)+1; e(3) = s(j1)+2;   % position in the global matrix
 e(4) = s(j2); e(5) = s(j2)+1; e(6) = s(j2)+2;
 e(7) = s(j3); e(8) = s(j3)+1; e(9) = s(j3)+2;
 e(10)= s(j4); e(11)= s(j5);   e(12)= s(j6);

 for j=1:12
   dof(j) = sol(e(j));
 end

 x1=p(j1,1); y1=p(j1,2);   % three vertices
 x2=p(j2,1); y2=p(j2,2);
 x3=p(j3,1); y3=p(j3,2);

 [a] = HCT_sys (x1,y1, x2,y2, x3,y3, dof);

 x=p(j4,1); y=p(j4,2);
 f(j4) = a(1)+a(2)*x+a(3)*y+a(4)*x^2+a(5)*x*y+a(6)*y^2 ...
        +a(7)*x^3+a(8)*x^2*y+a(9)*x*y^2+a(10)*y^3;
 x=p(j5,1); y=p(j5,2);
 f(j5) = a(11)+a(12)*x+a(13)*y+a(14)*x^2+a(15)*x*y+a(16)*y^2 ...
        +a(17)*x^3+a(18)*x^2*y+a(19)*x*y^2+a(20)*y^3;
 x=p(j6,1); y=p(j6,2);
 f(j6) = a(21)+a(22)*x+a(23)*y+a(24)*x^2+a(25)*x*y+a(26)*y^2 ...
        +a(27)*x^3+a(28)*x^2*y+a(29)*x*y^2+a(30)*y^3;
end

%for i=1:ne
%  f(c(i,4))=0.5*( f(c(i,1))+f(c(i,2)) );
%  f(c(i,5))=0.5*( f(c(i,2))+f(c(i,3)) );
%  f(c(i,6))=0.5*( f(c(i,3))+f(c(i,1)) );
%end

%-----
% extended connectivity matrix
% for 3-node sub-triangles
%-----

Ic=0;

for i=1:ne
 Ic=Ic+1;
 c3(Ic,1) = c(i,1); c3(Ic,2) = c(i,4); c3(Ic,3) = c(i,6);
 Ic=Ic+1;
 c3(Ic,1) = c(i,4); c3(Ic,2) = c(i,2); c3(Ic,3) = c(i,5);
 Ic=Ic+1;
 c3(Ic,1) = c(i,5); c3(Ic,2) = c(i,3); c3(Ic,3) = c(i,6);
 Ic=Ic+1;
 c3(Ic,1) = c(i,4); c3(Ic,2) = c(i,5); c3(Ic,3) = c(i,6);
end
                                                                                
trimesh (c3,p(:,1),p(:,2),f);
plot_6 (ne,ng,p,c,f);

%-----
% done
%-----

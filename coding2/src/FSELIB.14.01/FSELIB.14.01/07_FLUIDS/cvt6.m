clear all
close all

%===========
% Code cvt6
%
% Stokes flow in a cavity computed
% with 6-node triangular elements
%========================

%-----------
% input data
%-----------

visc = 1.0; % viscosity
V = 1.0;    % lid velocity 
NQ = 6;     % gauss-triangle quadrature
ndiv = 2;   % discretization level

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);

disp('Number of elements:'); ne

%-------------------------
% count the interior nodes
%-------------------------

inodes = 0;

for j=1:ng
 if(gfl(j,1)==0) 
  inodes = inodes+1;
 end
end

disp('Number of interior nodes:'); inodes

%------------------------------
% specify the boundary velocity
%------------------------------

for i=1:ng

 if(gfl(i,1)==1)

   gfl(i,2) = 0.0; % x velocity
   gfl(i,3) = 0.0; % y velocity

   if(p(i,2) > 0.99)
    gfl(i,2) = V;  % x velocity on the lid
   end

 end

end

%----------------------
% deform to a rectangle
%----------------------

def = 0.0; % example
def = 2/3; % example
def = 0.2; % example

for i=1:ng
  p(i,1) = p(i,1)*(1.0+def);
  p(i,2) = p(i,2)*(1.0-def);
end

%-------------------------------------
% assemble the global diffusion matrix
%          the Dx and Dy matrices
%-------------------------------------

gdm = zeros(ng,ng); % initialize
gDx = zeros(ne,ng); % initialize
gDy = zeros(ne,ng); % initialize

for l=1:ne          % loop over the elements

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);
j=c(l,4); x4=p(j,1); y4=p(j,2);
j=c(l,5); x5=p(j,1); y5=p(j,2);
j=c(l,6); x6=p(j,1); y6=p(j,2);

[edm_elm, arel] = edm6 ...
...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ);

[Dx, Dy] = cvt6_D ...
  ...
   (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ);

   for i=1:6
     i1 = c(l,i);
     for j=1:6
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
     gDx(l,i1) = Dx(i);
     gDy(l,i1) = Dy(i);
   end

end

%-------------------------
% compile the grand matrix
%-------------------------

for i=1:ng     % first block

 ngi = ng+i;

 for j=1:ng
  ngj = ng+j;
  Gm(i,j)   = visc*gdm(i,j);  Gm(i,ngj)   = 0.0;
  Gm(ngi,j) = 0.0;            Gm(ngi,ngj) = Gm(i,j);
 end

 for j=1:ne
  Gm(i,  2*ng+j) = -gDx(j,i);
  Gm(ngi,2*ng+j) = -gDy(j,i);
 end

end

% second block:

for i=1:ne

 for j=1:ng
  Gm(2*ng+i,j)    = -gDx(i,j);
  Gm(2*ng+i,ng+j) = -gDy(i,j);
 end
 for j=1:ne
  Gm(2*ng+i,2*ng+j) = 0;
 end

end

%------------
% system size
%------------

nsys = 2*ng+ne;

% disp (Gm);

%------------------------
% compute the right-hand side
% of the linear system
% and implement the Dirichlet BC
%------------------------

for i=1:nsys
 b(i) = 0.0;
end

for j=1:ng

 if(gfl(j,1)==1) 

   for i=1:nsys
    b(i) = b(i) - Gm(i,j)*gfl(j,2) - Gm(i,ng+j)*gfl(j,3);
    Gm(i,j) = 0; Gm(i,ng+j) = 0;
    Gm(j,i) = 0; Gm(ng+j,i) = 0;
   end

   Gm(j,j) = 1.0;
   Gm(ng+j,ng+j) = 1.0;
   b(j) = gfl(j,2);
   b(ng+j) = gfl(j,3);

 end

end

%disp (Gm);

%------------------------
% solve the linear system
%------------------------

Gm(:,nsys) = [];  % remove the last (nsys) column
Gm(nsys,:) = [];  % remove the last (nsys) row
 b(:,nsys) = [];  % remove the last (nsys) element

f = b/Gm'; 

%---------------------
% recover the velocity
%---------------------

for i=1:ng
 ux(i) = f(i);
 uy(i) = f(ng+i);
end

%---------------------
% recover the pressure
%---------------------

for i=1:ne-1
 press(i) = f(2*ng+i);
end
press(ne) =0.0;

%---------------------------------
% evaluate the nodal values of the
% by averaging pressure
%---------------------------------

cen = cen6(ne,ng,c);

for i=1:ng
   pg(i)=0.0;
   for j=2:cen(i,1)+1
    pg(i) = pg(i)+press(cen(i,j));
   end
   pg(i) = pg(i)/cen(i,1);
end

%-------------------------------
% plot the velocity vector field
%-------------------------------

figure(1)
hold on;
xlabel('x','fontsize',10)
ylabel('y','fontsize',10)
set(gca,'fontsize',15)
axis equal
quiver (p(:,1)',p(:,2)',ux,uy,'k');
%plot_6 (ne,ng,p,c,ux);
box on
axis((1+def)*[-1.0 1.0 -1.0, 1.0])

%-----------------------
% plot the element edges
%-----------------------

iplot=0;

if(iplot==1)

for i=1:ne

 i1=c(i,1);i2=c(i,2);i3=c(i,3);
 i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xp(1)=p(i1,1); xp(2)=p(i4,1); xp(3)=p(i2,1); xp(4)=p(i5,1);
                xp(5)=p(i3,1); xp(6)=p(i6,1); xp(7)=p(i1,1);
 yp(1)=p(i1,2); yp(2)=p(i4,2); yp(3)=p(i2,2); yp(4)=p(i5,2);
                yp(5)=p(i3,2); yp(6)=p(i6,2); yp(7)=p(i1,2);
 plot(xp, yp, 'k:');
%plot(xp, yp,'ko','markersize',5);

end

end

%---------------------------
% plot the boundary elements
%---------------------------

iplot=0;
iplot=1;

if(iplot==1)
                                                                                
for i=1:ne

 i1=c(i,1);i2=c(i,2);i3=c(i,3);
 i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xpb(1)=p(i1,1); xpb(2)=p(i4,1);
 ypb(1)=p(i1,2); ypb(2)=p(i4,2);
   if(efl(i,1)==1 & efl(i,4)==1) plot(xpb,ypb); end

 xpb(1)=p(i4,1); xpb(2)=p(i2,1);
 ypb(1)=p(i4,2); ypb(2)=p(i2,2);
   if(efl(i,4)==1 & efl(i,2)==1) plot(xpb,ypb); end

 xpb(1)=p(i2,1); xpb(2)=p(i5,1);
 ypb(1)=p(i2,2); ypb(2)=p(i5,2);
   if(efl(i,2)==1 & efl(i,5)==1) plot(xpb,ypb); end

 xpb(1)=p(i5,1); xpb(2)=p(i3,1);
 ypb(1)=p(i5,2); ypb(2)=p(i3,2);
   if(efl(i,5)==1 & efl(i,3)==1) plot(xpb,ypb); end
                                                                                
 xpb(1)=p(i3,1); xpb(2)=p(i6,1);
 ypb(1)=p(i3,2); ypb(2)=p(i6,2);
   if(efl(i,3)==1 & efl(i,6)==1) plot(xpb,ypb); end

 xpb(1)=p(i6,1); xpb(2)=p(i1,1);
 ypb(1)=p(i6,2); ypb(2)=p(i1,2);
   if(efl(i,6)==1 & efl(i,1)==1) plot(xpb,ypb); end

end

end

%--------------------------------
% plot the element pressure field
%--------------------------------

figure(2);
hold on;
xlabel('x','fontsize',10)
ylabel('y','fontsize',10)
zlabel('Pressure','fontsize',10)
set(gca,'fontsize',15)
%axis('equal')
box on
view(-30,26)

for i=1:ne

 i1=c(i,1);i2=c(i,2);i3=c(i,3);

 xplt(1) = (p(i1,1)+p(i2,1)+p(i3,1))/3.0; % element centroid
 yplt(1) = (p(i1,2)+p(i2,2)+p(i3,2))/3.0;
 zplt(1) = press(i);
 xplt(2) = xplt(1); yplt(2) = yplt(1);
 zplt(2) = 0.0;

 plot3(xplt,yplt,zplt,'k-','markersize',10)

end

%-----------------------
% plot the element edges
%-----------------------
                                                                                
for i=1:ne

 i1=c(i,1);i2=c(i,2);i3=c(i,3);
 i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xp(1)=p(i1,1); xp(2)=p(i4,1); xp(3)=p(i2,1); xp(4)=p(i5,1);
                xp(5)=p(i3,1); xp(6)=p(i6,1); xp(7)=p(i1,1);
 yp(1)=p(i1,2); yp(2)=p(i4,2); yp(3)=p(i2,2); yp(4)=p(i5,2);
                yp(5)=p(i3,2); yp(6)=p(i6,2); yp(7)=p(i1,2);
 plot(xp, yp, 'k:');
% plot(xp, yp,'o','markersize',10);
end

%------------------------------
% plot the nodal pressure field
%------------------------------

figure(3);
hold on
plot_6 (ne,ng,p,c,pg);
xlabel('x','fontsize',10)
ylabel('y','fontsize',10)
zlabel('Pressure','fontsize',10)
set(gca,'fontsize',15)
axis('equal')

%-----
% done
%-----

%--------------------------->
break

%------------------------
% plot the pressure field
%------------------------

figure; hold on;

for i=1:ne

 i1=c(i,1);i2=c(i,2);i3=c(i,3);
 i4=c(i,4);i5=c(i,5);i6=c(i,6);

 xelm(1) = p(i1,1);xelm(2) = p(i4,1);xelm(3) = p(i2,1);
 xelm(4) = p(i5,1);xelm(5) = p(i3,1);xelm(6) = p(i6,1);
 xelm(7) = p(i1,1);

 yelm(1) = p(i1,2);yelm(2) = p(i4,2);yelm(3) = p(i2,2);
 yelm(4) = p(i5,2);yelm(5) = p(i3,2);yelm(6) = p(i6,2);
 yelm(7) = p(i1,2);

 zelm(1) = press(i);zelm(2) = press(i);zelm(3) = press(i);
 zelm(4) = press(i);zelm(5) = press(i);zelm(6) = press(i);
 zelm(7) = press(i);

 plot3(xelm,yelm,zelm);

end

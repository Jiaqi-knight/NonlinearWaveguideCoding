
clear all
close all

%=========================================
% edges6
%
% visualize the edges of a 6-node triangle
%=========================================

figure(1)
hold on
axis equal
axis([-0.2 1 -0.5 1.6])
axis off

%----
for ipass=1:2  % one or two elements
%----

%-------------
% one triangle
%-------------

if(ipass==1)
 x1=-0.1; y1=-0.2; x2=0.93; y2= 0.1; x3=0.05; y3= 0.95; 
 x4= 0.8; y4=-0.1; x5=0.7;  y5= 0.2; x6=-0.1; y6=0.1; 
end

%-----------------
% another triangle
%-----------------

if(ipass==2)
 x1=0.7; y1= 1.0; x2=0.05; y2= 0.95; x3=0.93; y3= 0.1; 
 x4=0.5; y4= 1.0; x5=0.7;  y5= 0.2;  x6=0.90; y6= 0.5; 
end

%-----------------------
% relabel the nodes once
%-----------------------

% save
% x1s=x1; y1s=y1; x2s=x2; y2s=y2; x3s=x3; y3s=y3;
% x4s=x4; y4s=y4; x5s=x5; y5s=y5; x6s=x6; y6s=y6;
% relabel
% x1=x2s; y1=y2s; x2=x3s; y2=y3s; x3=x1s; y3=y1s;
% x4=x5s; y4=y5s; x5=x6s; y5=y6s; x6=x4s; y6=y4s;

%------------------------
% relabel the nodes twice
%------------------------

% save
% x1s=x1; y1s=y1; x2s=x2; y2s=y2; x3s=x3; y3s=y3;
% x4s=x4; y4s=y4; x5s=x5; y5s=y5; x6s=x6; y6s=y6;
% relabel
% x1=x2s; y1=y2s; x2=x3s; y2=y3s; x3=x1s; y3=y1s;
% x4=x5s; y4=y5s; x5=x6s; y5=y6s; x6=x4s; y6=y4s;

%---
for itest=1:2
%---

%---------------------------------------------
% compute the coefficients: alpha, beta, gamma
%---------------------------------------------

if(itest==1)

 [al, be, ga ] = elm6_abc ...
 ...
    (x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6);

elseif(itest==2)

 al=0.5; be=0.5; ga=0.5;

end

%--------
% prepare
%--------

np = 32;
step = 1.0/np;

alc = 1.0-al;
bec = 1.0-be;
gac = 1.0-ga;

%---------------
% plot the nodes
%---------------

plot (x1,y1,'ko','Markersize',10);
plot (x2,y2,'ko','Markersize',10);
plot (x3,y3,'ko','Markersize',10);
plot (x4,y4,'ko','Markersize',10);
plot (x5,y5,'ko','Markersize',10);
plot (x6,y6,'ko','Markersize',10);

%-----------
% edge 1-4-2
%-----------

for i=1:np+1
 xi=(i-1.0)*step;
 psi2 = xi*(xi-al)/alc ;
 psi4 = xi*(1.0-xi)/(al*alc) ;
 psi1 = 1.0-psi2-psi4;
 x(i) = x1*psi1 + x2*psi2 + x4*psi4;
 y(i) = y1*psi1 + y2*psi2 + y4*psi4;
end

if(itest==1)
 plot (x,y,'k');
elseif(itest==2)
 plot (x,y,'k:');
end

%-----------
% edge 1-6-3
%-----------

for i=1:np+1
 eta=(i-1.0)*step;
 psi3 = eta*(eta-be)/bec ;
 psi6 = eta*(1.0-eta)/(be*bec) ;
 psi1 = 1.0-psi3-psi6;
 x(i) = x1*psi1 + x3*psi3 + x6*psi6;
 y(i) = y1*psi1 + y3*psi3 + y6*psi6;
end

if(itest==1)
 plot (x,y,'k');
elseif(itest==2)
 plot (x,y,'k:');
end

%-----------
% edge 2-5-3
%-----------

for i=1:np+1
 xi=(i-1.0)*step;
 eta = 1.0-xi;
 psi2 = xi*(xi-al+eta*(al-ga)/gac)/alc;
 psi5 = xi*eta/(ga*gac);
 psi3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec;
 x(i) = x2*psi2 + x3*psi3 + x5*psi5;
 y(i) = y2*psi2 + y3*psi3 + y5*psi5;
end

if(itest==1)
 plot (x,y,'k');
elseif(itest==2)
 plot (x,y,'k:');
end

%-----
 end % of itest
%-----
 end % of ipass
%-----

%-----
% done
%-----


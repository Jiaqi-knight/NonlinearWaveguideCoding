
clear all
close all

%=========================================
% Driver for triangulating the unit sphere
% into three-node quadratic triangles
% by subdividing a regular icosahedron
%=========================================

 ndiv = 3;

%-------------
% triangulate
%-------------

 [npts,nelm,p,ne,n,nbe] = trgl3_octa (ndiv);

%---
% plot
%---

 figure(1)
 hold on
 axis equal
 box on
 axis off
 view(18,18)

 for i=1:nelm
  j1= n(i,1);
  j2= n(i,2);
  j3= n(i,3);
  px1 = p(j1,1); px2 = p(j2,1); px3 = p(j3,1);
  py1 = p(j1,2); py2 = p(j2,2); py3 = p(j3,2);
  pz1 = p(j1,3); pz2 = p(j2,3); pz3 = p(j3,3);
  patch([px1,px2,px3],[py1,py2,py3],[pz1,pz2,pz3],'w' ...
       ,'EdgeColor','k')
 end

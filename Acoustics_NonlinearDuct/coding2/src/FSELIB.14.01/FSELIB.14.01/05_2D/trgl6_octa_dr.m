close all
clear all

%==============
% generate a display a network of 6-node
% triangles descending from the icosahedron
%==============

ndiv = 2;

[npts,nelm,p,ne,n,nbe] = trgl6_octa (ndiv);

%----
% plot elements subdivided into 4 triangles
%----

figure(1)
hold on
axis equal
box on
axis off

%plot3(p(:,1),p(:,2),p(:,3),'ko')

for i=1:nelm
 j1= n(i,1);
 j2= n(i,2);
 j3= n(i,3);
 j4= n(i,4);
 j5= n(i,5);
 j6= n(i,6);
 px1 = p(j1,1); px2 = p(j2,1); px3 = p(j3,1);
 px4 = p(j4,1); px5 = p(j5,1); px6 = p(j6,1); 
 py1 = p(j1,2); py2 = p(j2,2); py3 = p(j3,2);
 py4 = p(j4,2); py5 = p(j5,2); py6 = p(j6,2); 
 pz1 = p(j1,3); pz2 = p(j2,3); pz3 = p(j3,3);
 pz4 = p(j4,3); pz5 = p(j5,3); pz6 = p(j6,3); 
%patch([px1,px4,px6],[py1,py4,py6],[pz1,pz4,pz6],'w')
%patch([px2,px5,px4],[py2,py5,py4],[pz2,pz5,pz4],'w')
%patch([px3,px6,px5],[py3,py6,py5],[pz3,pz6,pz5],'w')
%patch([px4,px5,px6],[py4,py5,py6],[pz4,pz5,pz6],'w')
 patch([px1,px4,px2,px5,px3,px6,px1], ...
       [py1,py4,py2,py5,py3,py6,py1], ...
       [pz1,pz4,pz2,pz5,pz3,pz6,pz1], 'w');
end

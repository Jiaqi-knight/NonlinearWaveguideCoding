
clear all
close all

%=================
% Script see_elm3
%
% scan the global nodes and display
% the host neighboring elements
%=================

ndiv=1;

%------
% prepare
%------

figure(1)
xlabel('x');
ylabel('y');
axis([-1.1 1.1 -1.1 1.1])

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);

%---------------------------------
% element-node connectivity matrix
%---------------------------------

cen = cen3(ne,ng,c);

for node=1:ng  % run over nodes

 node
 disp(gfl(node))
%clf

 for i=2:cen(node,1)+1  % run over host elements

   l = cen(node,i);    % element number

   j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2);
   j=c(l,2); xp(2)=p(j,1); yp(2)=p(j,2);
   j=c(l,3); xp(3)=p(j,1); yp(3)=p(j,2);
   j=c(l,1); xp(4)=p(j,1); yp(4)=p(j,2);


   cp(1)=1.0;  % another color code
   cp(2)=0.0;
   cp(3)=0.0;
   cp(4)=1.0;

   cp(1:4)=node; % color code

   patch(xp,yp,cp);
   hold on

 end
 plot(p(node,1),p(node,2),'ko')
 pause
 hold off
 plot(p(node,1),p(node,2),'ko')
end

%-----
% done
%-----

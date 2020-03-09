function [ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv)

%=================================
% triangulation of the unit disk
% into 3-node triangles
% by successive subdivisions
% of a 4-element parental structure
%
% LEGEND:
%
% efl: element flag
%      0 for an interior node
%      1 for a boundary node
% gfl: global node flag
%      0 for an interior node
%      1 for a boundary node
%==================================

%-----------------------------------
% parental structure with 4 elements
%-----------------------------------

iopt = 2;
iopt = 1;

if(iopt==1)

ne = 4;

x(1,1)= 0.0; y(1,1)= 0.0; efl(1,1)=0;  % first element
x(1,2)= 1.0; y(1,2)= 0.0; efl(1,2)=1;
x(1,3)= 0.0; y(1,3)= 1.0; efl(1,3)=1;

x(2,1)= 0.0; y(2,1)= 0.0; efl(2,1)=0;  % second element
x(2,2)= 0.0; y(2,2)= 1.0; efl(2,2)=1;
x(2,3)=-1.0; y(2,3)= 0.0; efl(2,3)=1;

x(3,1)= 0.0; y(3,1)= 0.0; efl(3,1)=0;  % third element
x(3,2)=-1.0; y(3,2)= 0.0; efl(3,2)=1;
x(3,3)= 0.0; y(3,3)=-1.0; efl(3,3)=1;

x(4,1)= 0.0; y(4,1)= 0.0; efl(4,1)=0;  % fourth element
x(4,2)= 0.0; y(4,2)=-1.0; efl(4,2)=1;
x(4,3)= 1.0; y(4,3)= 0.0; efl(4,3)=1;

elseif(iopt==2)

ne = 3;

bea=sin(2*pi/ne);

x(1,1)= 0.0; y(1,1)= 0.0; efl(1,1)=0;  % first element
x(1,2)= 1.0; y(1,2)= 0.0; efl(1,2)=1;
x(1,3)=-0.5; y(1,3)= bea; efl(1,3)=1;

x(2,1)= 0.0; y(2,1)= 0.0; efl(2,1)=0;  % second element
x(2,2)=-0.5; y(2,2)= bea; efl(2,2)=1;
x(2,3)=-0.5; y(2,3)=-bea; efl(2,3)=1;

x(3,1)= 0.0; y(3,1)= 0.0; efl(3,1)=0;  % third element
x(3,2)=-0.5; y(3,2)=-bea; efl(3,2)=1;
x(3,3)= 1.0; y(3,3)= 0.0; efl(3,3)=1;

end

%----------------
% refinement loop
%----------------

if(ndiv > 0)

for i=1:ndiv

 nm = 0; % count the new elements arising by each refinement loop
         % four elements are generated in each pass

 for j=1:ne   % loop over current elements

  % edge mid-points become vertex nodes

   x(j,4) = 0.5*(x(j,1)+x(j,2)); y(j,4) = 0.5*(y(j,1)+y(j,2));
   x(j,5) = 0.5*(x(j,2)+x(j,3)); y(j,5) = 0.5*(y(j,2)+y(j,3));
   x(j,6) = 0.5*(x(j,3)+x(j,1)); y(j,6) = 0.5*(y(j,3)+y(j,1));

   efl(j,4) = 0; if(efl(j,1)==1 & efl(j,2)==1) efl(j,4) = 1; end
   efl(j,5) = 0; if(efl(j,2)==1 & efl(j,3)==1) efl(j,5) = 1; end
   efl(j,6) = 0; if(efl(j,3)==1 & efl(j,1)==1) efl(j,6) = 1; end

  % assign vertex nodes to sub-elements
  % these will become the "new" elements

   nm = nm+1;
   xn(nm,1)=x(j,1); yn(nm,1)=y(j,1); efln(nm,1)=efl(j,1); %  first sub-element
   xn(nm,2)=x(j,4); yn(nm,2)=y(j,4); efln(nm,2)=efl(j,4);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  second sub-element
   xn(nm,2)=x(j,2); yn(nm,2)=y(j,2); efln(nm,2)=efl(j,2);
   xn(nm,3)=x(j,5); yn(nm,3)=y(j,5); efln(nm,3)=efl(j,5);
   nm = nm+1;
   xn(nm,1)=x(j,6); yn(nm,1)=y(j,6); efln(nm,1)=efl(j,6); %  third sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,3); yn(nm,3)=y(j,3); efln(nm,3)=efl(j,3);
   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  fourth sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);

 end % end of loop over current elements

%---

 ne = 4*ne;  % number of elements has increased
             % by a factor of 4

 for k=1:ne     % relabel the new points
                % and put them in the master list
                % project boundary nodes onto the unit circle
    for l=1:3
       x(k,l) =   xn(k,l);
       y(k,l) =   yn(k,l);
     efl(k,l) = efln(k,l);
      if(efl(k,l) == 1)
        rad = sqrt(x(k,l)^2+y(k,l)^2);
        x(k,l) = x(k,l)/rad;
        y(k,l) = y(k,l)/rad;
      end
    end

 end

end % end do of refinement loop
end % end if of refinement loop

%---------
% plotting
%---------

iplot = 1;
iplot = 0;

if(iplot==1)

 figure(88)
 hold on
 axis square
 xlabel('x','fontsize',15);
 ylabel('y','fontsize',15);
 set(gca,'fontsize',15)

 for i=1:ne
  xp(1)=x(i,1); xp(2)=x(i,2); xp(3)=x(i,3); xp(4)=x(i,1);
  yp(1)=y(i,1); yp(2)=y(i,2); yp(3)=y(i,3); yp(4)=y(i,1);
  plot(xp, yp,'-ko');
 end

end

%---------------------------------------------------
% define the global nodes and the connectivity table
%---------------------------------------------------

% three nodes of the first element are entered mannualy

p(1,1)=x(1,1); p(1,2)=y(1,1); gfl(1)=efl(1,1);
p(2,1)=x(1,2); p(2,2)=y(1,2); gfl(2)=efl(1,2);
p(3,1)=x(1,3); p(3,2)=y(1,3); gfl(3)=efl(1,3);

c(1,1) = 1;  % first  node of first element is global node 1
c(1,2) = 2;  % second node of first element is global node 2
c(1,3) = 3;  % third  node of first element is global node 3

ng = 3;

%---
% loop over further elements
% Iflag=0 will signal a new global node
%---

eps = 0.000001;  % tolerance

for i=2:ne        % loop over elements
 for j=1:3          % loop over element nodes

 Iflag = 0;

 for k=1:ng

  if(abs(x(i,j)-p(k,1)) < eps)
   if(abs(y(i,j)-p(k,2)) < eps)

     Iflag = 1;    % the node has been recorded previously
     c(i,j) = k;   % the jth local node of element i is the kth global node

   end
  end

 end

 if(Iflag==0)  % record the node
   ng = ng+1;
     p(ng,1) =   x(i,j);
     p(ng,2) =   y(i,j);
     gfl(ng) = efl(i,j);
   c(i,j) = ng;   % the jth local node of element is the new global node
 end

 end
end  % end of loop over elements

%-----
% done
%-----

return;

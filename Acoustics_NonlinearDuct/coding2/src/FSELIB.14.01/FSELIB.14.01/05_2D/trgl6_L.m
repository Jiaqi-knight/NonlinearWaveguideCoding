function [ne,ng,p,c,efl,gfl] = trgl6_L(ndiv)

%=========================================================
% Triangulation of an L-shaped domain into 6-node elements
% by the successive subdivision of 6 parental elements.
%========================================================

%-----------------------------------
% parental structure with 6 elements
%-----------------------------------

ne = 6;

x(1,1) = 0.0; y(1,1) = 0.0; efl(1,1) = 1;  % first element
x(1,2) = 1.0; y(1,2) = 0.0; efl(1,2) = 1;
x(1,3) = 0.0; y(1,3) = 1.0; efl(1,3) = 1;

x(1,4) = 0.5*(x(1,1)+x(1,2));
y(1,4) = 0.5*(y(1,1)+y(1,2)); efl(1,4) = 0;
x(1,5) = 0.5*(x(1,2)+x(1,3));
y(1,5) = 0.5*(y(1,2)+y(1,3)); efl(1,5) = 0;
x(1,6) = 0.5*(x(1,3)+x(1,1));
y(1,6) = 0.5*(y(1,3)+y(1,1)); efl(1,6) = 1;

x(2,1) = 1.0; y(2,1) = 0.0; efl(2,1) = 1;  % second element
x(2,2) = 1.0; y(2,2) = 1.0; efl(2,2) = 1;
x(2,3) = 0.0; y(2,3) = 1.0; efl(2,3) = 1;

x(2,4) = 0.5*(x(2,1)+x(2,2));
y(2,4) = 0.5*(y(2,1)+y(2,2)); efl(2,4) = 1;
x(2,5) = 0.5*(x(2,2)+x(2,3));
y(2,5) = 0.5*(y(2,2)+y(2,3)); efl(2,5) = 1;
x(2,6) = 0.5*(x(2,3)+x(2,1));
y(2,6) = 0.5*(y(2,3)+y(2,1)); efl(2,6) = 0;

x(3,1) = 0.0; y(3,1) = 0.0; efl(3,1) = 1;  % third element
x(3,2) =-1.0; y(3,2) = 0.0; efl(3,2) = 1;
x(3,3) = 0.0; y(3,3) =-1.0; efl(3,3) = 1;

x(3,4) = 0.5*(x(3,1)+x(3,2));
y(3,4) = 0.5*(y(3,1)+y(3,2)); efl(3,4) = 1;
x(3,5) = 0.5*(x(3,2)+x(3,3));
y(3,5) = 0.5*(y(3,2)+y(3,3)); efl(3,5) = 0;
x(3,6) = 0.5*(x(3,3)+x(3,1));
y(3,6) = 0.5*(y(3,3)+y(3,1)); efl(3,6) = 0;

x(4,1) =-1.0; y(4,1) = 0.0; efl(4,1)=1;  % fourth element
x(4,2) =-1.0; y(4,2) =-1.0; efl(4,2)=1;
x(4,3) = 0.0; y(4,3) =-1.0; efl(4,3)=1;

x(4,4) = 0.5*(x(4,1)+x(4,2));
y(4,4) = 0.5*(y(4,1)+y(4,2));efl(4,4)=1;
x(4,5) = 0.5*(x(4,2)+x(4,3));
y(4,5) = 0.5*(y(4,2)+y(4,3));efl(4,5)=1;
x(4,6) = 0.5*(x(4,3)+x(4,1));
y(4,6) = 0.5*(y(4,3)+y(4,1));efl(4,6)=0;

x(5,1) = 0.0; y(5,1) = 0.0; efl(5,1)=0;  % fifth element
x(5,2) = 0.0; y(5,2) =-1.0; efl(5,2)=1;
x(5,3) = 1.0; y(5,3) = 0.0; efl(5,3)=1;

x(5,4) = 0.5*(x(5,1)+x(5,2));
y(5,4) = 0.5*(y(5,1)+y(5,2));efl(5,4)=0;
x(5,5) = 0.5*(x(5,2)+x(5,3));
y(5,5) = 0.5*(y(5,2)+y(5,3));efl(5,5)=0;
x(5,6) = 0.5*(x(5,3)+x(5,1));
y(5,6) = 0.5*(y(5,3)+y(5,1));efl(5,6)=0;

x(6,1) = 1.0; y(6,1) = 0.0; efl(6,1)=1;  % sixth element
x(6,2) = 0.0; y(6,2) =-1.0; efl(6,2)=1;
x(6,3) = 1.0; y(6,3) =-1.0; efl(6,3)=1;
x(6,4) = 0.5*(x(6,1)+x(6,2));
y(6,4) = 0.5*(y(6,1)+y(6,2));efl(6,4)=0;
x(6,5) = 0.5*(x(6,2)+x(6,3));
y(6,5) = 0.5*(y(6,2)+y(6,3));efl(6,5)=1;
x(6,6) = 0.5*(x(6,3)+x(6,1));
y(6,6) = 0.5*(y(6,3)+y(6,1));efl(6,6)=1;

% ---------------
% refinement loop
% ---------------

if(ndiv > 0)

for i=1:ndiv

 nm = 0; % count the new elements arising by each refinement loop
         % four elements will be generated in each pass

 for j=1:ne   % loop over current elements

  % assign vertex nodes to sub-elements
  % these will become the "new" elements

   nm = nm+1;
   xn(nm,1)=x(j,1); yn(nm,1)=y(j,1); efln(nm,1)=efl(j,1); %  first sub-element
   xn(nm,2)=x(j,4); yn(nm,2)=y(j,4); efln(nm,2)=efl(j,4);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end

   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  second sub-element
   xn(nm,2)=x(j,2); yn(nm,2)=y(j,2); efln(nm,2)=efl(j,2);
   xn(nm,3)=x(j,5); yn(nm,3)=y(j,5); efln(nm,3)=efl(j,5);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end

   nm = nm+1;
   xn(nm,1)=x(j,6); yn(nm,1)=y(j,6); efln(nm,1)=efl(j,6); %  third sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,3); yn(nm,3)=y(j,3); efln(nm,3)=efl(j,3);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end

   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  fourth sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end

 end % end of loop over current elements

%---

 ne = 4*ne;  % number of elements has increased
             % by a factor of four

 for k=1:ne     % relabel the new points
                % and put them in the master list
    for l=1:6
       x(k,l)=xn(k,l);
       y(k,l)=yn(k,l);
     efl(k,l)=efln(k,l);
    end

 end

 for k=1:ne     % reposition the mid-nodes
    if(efl(k,4) == 0)
      x(k,4)=0.5*(x(k,1)+x(k,2)); y(k,4)=0.5*(y(k,1)+y(k,2)); end
    if(efl(k,5) == 0)
      x(k,5)=0.5*(x(k,2)+x(k,3)); y(k,5)=0.5*(y(k,2)+y(k,3)); end
    if(efl(k,6) == 0)
      x(k,6)=0.5*(x(k,3)+x(k,1)); y(k,6)=0.5*(y(k,3)+y(k,1)); end
 end

end % do of refinement loop
end % if of refinement loop

%-----
% reset the boundary flags
%-----

 for k=1:ne     % reposition the mid-nodes
   for l=1:6
    if(efl(k,l)==1)
     test= x(k,l)*(x(k,l)-1.0)*(x(k,l)+1.0) ...
          *y(k,l)*(y(k,l)-1.0)*(y(k,l)+1.0);
     if(abs(test)>0.001)
      efl(k,l) = 0;
     end
   end
   end
 end

% --------
% plotting
% --------

iplot = 1;
iplot = 0;

if(iplot==1)

 figure(88)
 hold on
 xlabel('x');
 ylabel('y');
 for i=1:ne
  xp(1)=x(i,1);xp(2)=x(i,4);xp(3)=x(i,2);xp(4)=x(i,5);
               xp(5)=x(i,3);xp(6)=x(i,6);xp(7)=x(i,1);
  yp(1)=y(i,1);yp(2)=y(i,4);yp(3)=y(i,2);
               yp(4)=y(i,5);yp(5)=y(i,3);yp(6)=y(i,6);yp(7)=y(i,1);
  plot(xp, yp,'-o');
 end

end

% --------
% define the global nodes and the connectivity table
% --------

% six nodes of the first element are entered mannualy

p(1,1)=x(1,1); p(1,2)=y(1,1); gfl(1)=efl(1,1);
p(2,1)=x(1,2); p(2,2)=y(1,2); gfl(2)=efl(1,2);
p(3,1)=x(1,3); p(3,2)=y(1,3); gfl(3)=efl(1,3);
p(4,1)=x(1,4); p(4,2)=y(1,4); gfl(4)=efl(1,4);
p(5,1)=x(1,5); p(5,2)=y(1,5); gfl(5)=efl(1,5);
p(6,1)=x(1,6); p(6,2)=y(1,6); gfl(6)=efl(1,6);

c(1,1) = 1;  % first  node of first element is global node 1
c(1,2) = 2;  % second node of first element is global node 2
c(1,3) = 3;  % third  node of first element is global node 3
c(1,4) = 4;  % fourth node of first element is global node 4
c(1,5) = 5;  % fifth  node of first element is global node 5
c(1,6) = 6;  % sixth  node of first element is global node 6

ng = 6;

% loop over further elements

% Iflag=0 will signal a new global node

eps = 0.000001;

for i=2:ne        % loop over elements
 for j=1:6          % loop over element nodes

 Iflag=0;
 for k=1:ng

  if(abs(x(i,j)-p(k,1)) < eps)
   if(abs(y(i,j)-p(k,2)) < eps)
     Iflag = 1;    % the node has been recorded previously
     c(i,j) = k;   % the jth local node of element i is the kth global node
   end
  end

 end

 if(Iflag==0)  % record a new node
   ng = ng+1;
   p(ng,1) = x(i,j);
   p(ng,2) = y(i,j);
   gfl(ng) = efl(i,j);
   c(i,j) = ng;   % the jth local node of element is the new global node
 end

 end
end  % end of loop over elements

%-----
% done
%-----

return;

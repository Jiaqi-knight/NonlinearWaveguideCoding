function [ne,ng,p,c,efl,gfl] = trgl6_rc (a,ndiv)

%=======================================================
% FSELIB
%
% Triangulation of a rectangle with a circular hole
% of radius ''a'' into 6-node elements
% by the successive subdivision of 20 parental elements.
%=======================================================

%---------
% parental structure with 20 elements
%---------

ne = 20;

x(1,1) =   a;y(1,1) =  -a; efl(1,1)=2;  % first element
x(1,2) = 1.0;y(1,2) =-1.0; efl(1,2)=1;
x(1,3) = 1.0;y(1,3) = 0.0; efl(1,3)=0;
x(1,4) = 0.5*(x(1,1)+x(1,2)); y(1,4) = 0.5*(y(1,1)+y(1,2));efl(1,4)=0;
x(1,5) = 0.5*(x(1,2)+x(1,3)); y(1,5) = 0.5*(y(1,2)+y(1,3));efl(1,5)=0;
x(1,6) = 0.5*(x(1,3)+x(1,1)); y(1,6) = 0.5*(y(1,3)+y(1,1));efl(1,6)=0;

x(2,1) =   a;y(2,1) =  -a; efl(2,1)=2;  % second element
x(2,2) = 1.0;y(2,2) = 0.0; efl(2,2)=0;
x(2,3) =   a;y(2,3) =   a; efl(2,3)=2;
x(2,4) = 0.5*(x(2,1)+x(2,2));y(2,4) = 0.5*(y(2,1)+y(2,2));efl(2,4)=0;
x(2,5) = 0.5*(x(2,2)+x(2,3));y(2,5) = 0.5*(y(2,2)+y(2,3));efl(2,5)=0;
x(2,6) = 0.5*(x(2,3)+x(2,1));y(2,6) = 0.5*(y(2,3)+y(2,1));efl(2,6)=2;

x(3,1) =   a;y(3,1) =   a; efl(3,1)=2;  % third element
x(3,2) = 1.0;y(3,2) = 0.0; efl(3,2)=0;
x(3,3) = 1.0;y(3,3) = 1.0; efl(3,3)=1;
x(3,4) = 0.5*(x(3,1)+x(3,2));y(3,4) = 0.5*(y(3,1)+y(3,2));efl(3,4)=0;
x(3,5) = 0.5*(x(3,2)+x(3,3));y(3,5) = 0.5*(y(3,2)+y(3,3));efl(3,5)=0;
x(3,6) = 0.5*(x(3,3)+x(3,1));y(3,6) = 0.5*(y(3,3)+y(3,1));efl(3,6)=0;

x(4,1) =   a;y(4,1) =   a; efl(4,1)=2;  % fourth element
x(4,2) = 1.0;y(4,2) = 1.0; efl(4,2)=1;
x(4,3) = 0.0;y(4,3) = 1.0; efl(4,3)=1;
x(4,4) = 0.5*(x(4,1)+x(4,2));y(4,4) = 0.5*(y(4,1)+y(4,2));efl(4,4)=0;
x(4,5) = 0.5*(x(4,2)+x(4,3));y(4,5) = 0.5*(y(4,2)+y(4,3));efl(4,5)=1;
x(4,6) = 0.5*(x(4,3)+x(4,1));y(4,6) = 0.5*(y(4,3)+y(4,1));efl(4,6)=0;

x(5,1) = 0.0;y(5,1) = 1.0; efl(5,1)=1;  % fifth element
x(5,2) =  -a;y(5,2) =   a; efl(5,2)=2;
x(5,3) =   a;y(5,3) =   a; efl(5,3)=2;
x(5,4) = 0.5*(x(5,1)+x(5,2));y(5,4) = 0.5*(y(5,1)+y(5,2));efl(5,4)=0;
x(5,5) = 0.5*(x(5,2)+x(5,3));y(5,5) = 0.5*(y(5,2)+y(5,3));efl(5,5)=2;
x(5,6) = 0.5*(x(5,3)+x(5,1));y(5,6) = 0.5*(y(5,3)+y(5,1));efl(5,6)=0;

x(6,1) =  -a;y(6,1) =   a; efl(6,1)=2;  % sixth element
x(6,2) = 0.0;y(6,2) = 1.0; efl(6,2)=1;
x(6,3) =-1.0;y(6,3) = 1.0; efl(6,3)=1;
x(6,4) = 0.5*(x(6,1)+x(6,2));y(6,4) = 0.5*(y(6,1)+y(6,2));efl(6,4)=0;
x(6,5) = 0.5*(x(6,2)+x(6,3));y(6,5) = 0.5*(y(6,2)+y(6,3));efl(6,5)=1;
x(6,6) = 0.5*(x(6,3)+x(6,1));y(6,6) = 0.5*(y(6,3)+y(6,1));efl(6,6)=0;

%---
% reflect to generate elements: 7-12
%---

for i=1:6
 for j=1:6
  x(6+i,j)=-x(i,j);y(6+i,j)=-y(i,j);efl(6+i,j)=efl(i,j);
 end
end

efl(7,1) =2;efl(7,2)=1;efl(7,3)=1;efl(7,4)=0;efl(7,5)=1;efl(7,6)=0;
efl(8,1) =2;efl(8,2)=1;efl(8,3)=2;efl(8,4)=0;efl(8,5)=0;efl(8,6)=2;
efl(9,1) =2;efl(9,2)=1;efl(9,3)=1;efl(9,4)=0;efl(9,5)=1;efl(9,6)=0;
efl(10,1)=2;efl(10,2)=1;efl(10,3)=1;efl(10,4)=0;efl(10,5)=1;efl(10,6)=0;
efl(11,1)=1;efl(11,2)=2;efl(11,3)=2;efl(11,4)=0;efl(11,5)=2;efl(11,6)=0;
efl(12,1)=2;efl(12,2)=1;efl(12,3)=1;efl(12,4)=0;efl(12,5)=1;efl(12,6)=0;

%---
% project on the circle
%---

for  i=1:12
  for j=1:6
    if (efl(i,j)==2)
     rad = a/sqrt(x(i,j)^2+y(i,j)^2);
     x(i,j) = x(i,j)*rad;y(i,j) = y(i,j)*rad;
    end
 end
end

%---
% elements on the right part
%---

x(13,1) = 1.0; y(13,1) =  -1.0; efl(13,1)=1;  % element 13
x(13,2) = 2.0; y(13,2) =  -1.0; efl(13,2)=1; 
x(13,3) = 1.0; y(13,3) =   0.0; efl(13,3)=0;
x(13,4) = 0.5*(x(13,1)+x(13,2));y(13,4) = 0.5*(y(13,1)+y(13,2));efl(13,4)=1;
x(13,5) = 0.5*(x(13,2)+x(13,3));y(13,5) = 0.5*(y(13,2)+y(13,3));efl(13,5)=0;
x(13,6) = 0.5*(x(13,3)+x(13,1));y(13,6) = 0.5*(y(13,3)+y(13,1));efl(13,6)=0;

x(14,1) = 2.0;y(14,1) =  -1.0; efl(14,1)=1;  % element 14
x(14,2) = 2.0;y(14,2) =  -0.0; efl(14,2)=0;
x(14,3) = 1.0;y(14,3) =   0.0; efl(14,3)=0;
x(14,4) = 0.5*(x(14,1)+x(14,2));y(14,4) = 0.5*(y(14,1)+y(14,2));efl(14,4)=0;
x(14,5) = 0.5*(x(14,2)+x(14,3));y(14,5) = 0.5*(y(14,2)+y(14,3));efl(14,5)=0;
x(14,6) = 0.5*(x(14,3)+x(14,1));y(14,6) = 0.5*(y(14,3)+y(14,1));efl(14,6)=0;

x(15,1) = 2.0;y(15,1) =   0.0; efl(15,1)=0;  % element 15
x(15,2) = 2.0;y(15,2) =   1.0; efl(15,2)=1;
x(15,3) = 1.0;y(15,3) =   0.0; efl(15,3)=0; 
x(15,4) = 0.5*(x(15,1)+x(15,2));y(15,4) = 0.5*(y(15,1)+y(15,2));efl(15,4)=0;
x(15,5) = 0.5*(x(15,2)+x(15,3));y(15,5) = 0.5*(y(15,2)+y(15,3));efl(15,5)=0;
x(15,6) = 0.5*(x(15,3)+x(15,1));y(15,6) = 0.5*(y(15,3)+y(15,1));efl(15,6)=0;

x(16,1) = 1.0;y(16,1) =   0.0; efl(16,1)=0;  % element 16
x(16,2) = 2.0;y(16,2) =   1.0; efl(16,2)=1; 
x(16,3) = 1.0;y(16,3) =   1.0; efl(16,3)=1; 
x(16,4) = 0.5*(x(16,1)+x(16,2));y(16,4) = 0.5*(y(16,1)+y(16,2));efl(16,4)=0;
x(16,5) = 0.5*(x(16,2)+x(16,3));y(16,5) = 0.5*(y(16,2)+y(16,3));efl(16,5)=1;
x(16,6) = 0.5*(x(16,3)+x(16,1));y(16,6) = 0.5*(y(16,3)+y(16,1));efl(16,6)=0;

x(17,1) = 2.0;y(17,1) =  -1.0; efl(17,1)=1;  % element 17
x(17,2) = 3.0;y(17,2) =  -1.0; efl(17,2)=1;
x(17,3) = 2.0;y(17,3) =   0.0; efl(17,3)=0;
x(17,4) = 0.5*(x(17,1)+x(17,2));y(17,4) = 0.5*(y(17,1)+y(17,2));efl(17,4)=1;
x(17,5) = 0.5*(x(17,2)+x(17,3));y(17,5) = 0.5*(y(17,2)+y(17,3));efl(17,5)=0;
x(17,6) = 0.5*(x(17,3)+x(17,1));y(17,6) = 0.5*(y(17,3)+y(17,1));efl(17,6)=0;

x(18,1) = 3.0;y(18,1) =  -1.0; efl(18,1)=1;  % element 18
x(18,2) = 3.0;y(18,2) =   0.0; efl(18,2)=1;
x(18,3) = 2.0;y(18,3) =   0.0; efl(18,3)=0; 
x(18,4) = 0.5*(x(18,1)+x(18,2));y(18,4) = 0.5*(y(18,1)+y(18,2));efl(18,4)=1;
x(18,5) = 0.5*(x(18,2)+x(18,3));y(18,5) = 0.5*(y(18,2)+y(18,3));efl(18,5)=0;
x(18,6) = 0.5*(x(18,3)+x(18,1));y(18,6) = 0.5*(y(18,3)+y(18,1));efl(18,6)=0;

x(19,1) = 3.0;y(19,1) =   0.0; efl(19,1)=1;  % element 19
x(19,2) = 3.0;y(19,2) =   1.0; efl(19,2)=1;
x(19,3) = 2.0;y(19,3) =   0.0; efl(19,3)=0;
x(19,4) = 0.5*(x(19,1)+x(19,2));y(19,4) = 0.5*(y(19,1)+y(19,2));efl(19,4)=1;
x(19,5) = 0.5*(x(19,2)+x(19,3));y(19,5) = 0.5*(y(19,2)+y(19,3));efl(19,5)=0;
x(19,6) = 0.5*(x(19,3)+x(19,1));y(19,6) = 0.5*(y(19,3)+y(19,1));efl(19,6)=0;

x(20,1) = 3.0;y(20,1) =   1.0; efl(20,1)=1;  % element 20
x(20,2) = 2.0;y(20,2) =   1.0; efl(20,2)=1; 
x(20,3) = 2.0;y(20,3) =   0.0; efl(20,3)=0; 
x(20,4) = 0.5*(x(20,1)+x(20,2));y(20,4) = 0.5*(y(20,1)+y(20,2));efl(20,4)=1;
x(20,5) = 0.5*(x(20,2)+x(20,3));y(20,5) = 0.5*(y(20,2)+y(20,3));efl(20,5)=0;
x(20,6) = 0.5*(x(20,3)+x(20,1));y(20,6) = 0.5*(y(20,3)+y(20,1));efl(20,6)=0;

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
                   if(efln(nm,1)==2 & efln(nm,2)==2) efln(nm,4) = 2; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
                   if(efln(nm,2)==2 & efln(nm,3)==2) efln(nm,5) = 2; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end
                   if(efln(nm,3)==2 & efln(nm,1)==2) efln(nm,6) = 2; end

   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  second sub-element
   xn(nm,2)=x(j,2); yn(nm,2)=y(j,2); efln(nm,2)=efl(j,2);
   xn(nm,3)=x(j,5); yn(nm,3)=y(j,5); efln(nm,3)=efl(j,5);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
                   if(efln(nm,1)==2 & efln(nm,2)==2) efln(nm,4) = 2; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
                   if(efln(nm,2)==2 & efln(nm,3)==2) efln(nm,5) = 2; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end
                   if(efln(nm,3)==2 & efln(nm,1)==2) efln(nm,6) = 2; end

   nm = nm+1;
   xn(nm,1)=x(j,6); yn(nm,1)=y(j,6); efln(nm,1)=efl(j,6); %  third sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,3); yn(nm,3)=y(j,3); efln(nm,3)=efl(j,3);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
                   if(efln(nm,1)==2 & efln(nm,2)==2) efln(nm,4) = 2; end
   efln(nm,5) = 0; if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
                   if(efln(nm,2)==2 & efln(nm,3)==2) efln(nm,5) = 2; end
   efln(nm,6) = 0; if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end
                   if(efln(nm,3)==2 & efln(nm,1)==2) efln(nm,6) = 2; end

   nm = nm+1;
   xn(nm,1)=x(j,4); yn(nm,1)=y(j,4); efln(nm,1)=efl(j,4); %  fourth sub-element
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,6); yn(nm,3)=y(j,6); efln(nm,3)=efl(j,6);
   xn(nm,4) = 0.5*(xn(nm,1)+xn(nm,2));yn(nm,4) = 0.5*(yn(nm,1)+yn(nm,2));
   xn(nm,5) = 0.5*(xn(nm,2)+xn(nm,3));yn(nm,5) = 0.5*(yn(nm,2)+yn(nm,3));
   xn(nm,6) = 0.5*(xn(nm,3)+xn(nm,1));yn(nm,6) = 0.5*(yn(nm,3)+yn(nm,1));

   efln(nm,4) = 0; efln(nm,5) = 0; efln(nm,6) = 0; 
   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,4) = 1; end
   if(efln(nm,1)==2 & efln(nm,2)==2) efln(nm,4) = 2; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,5) = 1; end
   if(efln(nm,2)==2 & efln(nm,3)==2) efln(nm,5) = 2; end
   if(efln(nm,3)==1 & efln(nm,1)==1) efln(nm,6) = 1; end
   if(efln(nm,3)==2 & efln(nm,1)==2) efln(nm,6) = 2; end

 end % end of loop over current elements

%---

 ne = 4*ne;  % number of elements has increased
             % by a factor of four

 for k=1:ne     % relabel the new points
                % and put them in the master list
    for l=1:6
     x(k,l)=xn(k,l); y(k,l)=yn(k,l);  efl(k,l)=efln(k,l);
     if (efl(k,l)==2)
      rad = a/sqrt(x(k,l)^2+y(k,l)^2);
      x(k,l) = x(k,l)*rad;y(k,l) = y(k,l)*rad;
     end
    end
 end

end % end do of refinement loop
end % end if of refinement loop

%----------
% plotting
%---------

for i=1:ne
 xp(1)=x(i,1);xp(2)=x(i,4);xp(3)=x(i,2);xp(4)=x(i,5);
              xp(5)=x(i,3);xp(6)=x(i,6);xp(7)=x(i,1);
 yp(1)=y(i,1);yp(2)=y(i,4);yp(3)=y(i,2);
              yp(4)=y(i,5);yp(5)=y(i,3);yp(6)=y(i,6);yp(7)=y(i,1);
%plot(xp, yp,'o-'); hold on
end
%axis([-1 3 -2 2]);
%xlabel('x');ylabel('y');

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

 if(Iflag==0)  % record the node
   ng = ng+1;
   p(ng,1) = x(i,j);
   p(ng,2) = y(i,j);
   gfl(ng) = efl(i,j);
   c(i,j) = ng;   % the jth local node of element is the new global node
 end

 end
end  % end of loop over elements

%-----
% Done
%-----

return;

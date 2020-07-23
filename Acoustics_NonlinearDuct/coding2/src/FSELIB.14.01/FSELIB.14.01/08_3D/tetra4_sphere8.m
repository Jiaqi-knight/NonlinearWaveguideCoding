function [ne,ng,p,c,efl,gfl] = tetra4_sphere (ndiv)

%=====================================================
% Discretize a sphere into 4-node tetrahedral elements
% by successive subdivisions of an 8-element
% parental structure
% 
% LEGEND:
%
% efl: element flag
% gfl: global node flag
%====================================

act = 1.0;

%---------------------------------------
% parental structure with eight elements
%---------------------------------------

ne = 8;

x(1,1)= 0.0; y(1,1)= 0.0; z(1,1)= 0.0; efl(1,1)=0;  % first element
x(1,2)= act; y(1,2)= 0.0; z(1,2)= 0.0; efl(1,2)=1;
x(1,3)= 0.0; y(1,3)= act; z(1,3)= 0.0; efl(1,3)=1;
x(1,4)= 0.0; y(1,4)= 0.0; z(1,4)= act; efl(1,4)=1;

x(2,1)= 0.0; y(2,1)= 0.0; z(2,1)= 0.0; efl(2,1)=0;  % second element
x(2,2)= 0.0; y(2,2)=-act; z(2,2)= 0.0; efl(2,2)=1;
x(2,3)= act; y(2,3)= 0.0; z(2,3)= 0.0; efl(2,3)=1;
x(2,4)= 0.0; y(2,4)= 0.0; z(2,4)= act; efl(2,4)=1;

x(3,1)= 0.0; y(3,1)= 0.0; z(3,1)= 0.0; efl(3,1)=0;  % third element
x(3,2)= 0.0; y(3,2)= act; z(3,2)= 0.0; efl(3,2)=1;
x(3,3)=-act; y(3,3)= 0.0; z(3,3)= 0.0; efl(3,3)=1;
x(3,4)= 0.0; y(3,4)= 0.0; z(3,4)= act; efl(3,4)=1;

x(4,1)= 0.0; y(4,1)= 0.0; z(4,1)= 0.0; efl(4,1)=0;  % fourth element
x(4,2)=-act; y(4,2)= 0.0; z(4,2)= 0.0; efl(4,2)=1;
x(4,3)= 0.0; y(4,3)=-act; z(4,3)= 0.0; efl(4,3)=1;
x(4,4)= 0.0; y(4,4)= 0.0; z(4,4)= act; efl(4,4)=1;

x(5,1)= 0.0; y(5,1)= 0.0; z(5,1)= 0.0; efl(5,1)=0;  % fifth element
x(5,2)= 0.0; y(5,2)= act; z(5,2)= 0.0; efl(5,3)=1;
x(5,3)= act; y(5,3)= 0.0; z(5,3)= 0.0; efl(5,2)=1;
x(5,4)= 0.0; y(5,4)= 0.0; z(5,4)=-act; efl(5,4)=1;

x(6,1)= 0.0; y(6,1)= 0.0; z(6,1)= 0.0; efl(6,1)=0;  % sixth element
x(6,2)= act; y(6,3)= 0.0; z(6,2)= 0.0; efl(6,2)=1;
x(6,3)= 0.0; y(6,3)=-act; z(6,3)= 0.0; efl(6,3)=1;
x(6,4)= 0.0; y(6,4)= 0.0; z(6,4)=-act; efl(6,4)=1;

x(7,1)= 0.0; y(7,1)= 0.0; z(7,1)= 0.0; efl(7,1)=0;  % seventh element
x(7,2)=-act; y(7,2)= 0.0; z(7,2)= 0.0; efl(7,2)=1;
x(7,3)= 0.0; y(7,3)= act; z(7,3)= 0.0; efl(7,3)=1;
x(7,4)= 0.0; y(7,4)= 0.0; z(7,4)=-act; efl(7,4)=1;

x(8,1)= 0.0; y(8,1)= 0.0; z(8,1)= 0.0; efl(8,1)=0;  % eighth element
x(8,2)= 0.0; y(8,2)=-act; z(8,2)= 0.0; efl(8,2)=1;
x(8,3)=-act; y(8,3)= 0.0; z(8,3)= 0.0; efl(8,3)=1;
x(8,4)= 0.0; y(8,4)= 0.0; z(8,4)=-act; efl(8,4)=1;

%----------------
% refinement loop
%----------------

if(ndiv > 0)

for i=1:ndiv

 nm = 0; % count the new elements arising by each refinement loop
         % eight elements will be generated in each pass

 for j=1:ne   % loop over current elements

   % edge mid-nodes

   x12 = 0.5*(x(j,1)+x(j,2));
   y12 = 0.5*(y(j,1)+y(j,2));
   z12 = 0.5*(z(j,1)+z(j,2));
   efl12 = 0; if(efl(j,1)==1 & efl(j,2)==1) efl12 = 1; end

   x13 = 0.5*(x(j,1)+x(j,3));
   y13 = 0.5*(y(j,1)+y(j,3));
   z13 = 0.5*(z(j,1)+z(j,3));
   efl13 = 0; if(efl(j,1)==1 & efl(j,3)==1) efl13 = 1; end

   x14 = 0.5*(x(j,1)+x(j,4));
   y14 = 0.5*(y(j,1)+y(j,4));
   z14 = 0.5*(z(j,1)+z(j,4));
   efl14 = 0; if(efl(j,1)==1 & efl(j,4)==1) efl14 = 1; end

   x23 = 0.5*(x(j,2)+x(j,3));
   y23 = 0.5*(y(j,2)+y(j,3));
   z23 = 0.5*(z(j,2)+z(j,3));
   efl23 = 0; if(efl(j,2)==1 & efl(j,3)==1) efl23 = 1; end

   x24 = 0.5*(x(j,2)+x(j,4));
   y24 = 0.5*(y(j,2)+y(j,4));
   z24 = 0.5*(z(j,2)+z(j,4));
   efl24 = 0; if(efl(j,2)==1 & efl(j,4)==1) efl24 = 1; end

   x34 = 0.5*(x(j,3)+x(j,4));
   y34 = 0.5*(y(j,3)+y(j,4));
   z34 = 0.5*(z(j,3)+z(j,4));
   efl34 = 0; if(efl(j,3)==1 & efl(j,4)==1) efl34 = 1; end

  % assign vertex nodes to sub-elements
  % these will become the "new" elements

   nm = nm+1;  % first subelement

   xn(nm,1)=x13;    yn(nm,1)=y13;    zn(nm,1)=z13;    efln(nm,1)=efl13;
   xn(nm,2)=x23;    yn(nm,2)=y23;    zn(nm,2)=z23;    efln(nm,2)=efl23;
   xn(nm,3)=x(j,3); yn(nm,3)=y(j,3); zn(nm,3)=z(j,3); efln(nm,3)=efl(j,3);
   xn(nm,4)=x34;    yn(nm,4)=y34;    zn(nm,4)=z34;    efln(nm,4)=efl34;

   nm = nm+1;  % second subelement

   xn(nm,1)=x12;    yn(nm,1)=y12;    zn(nm,1)=z12;    efln(nm,1)=efl12;
   xn(nm,2)=x(j,2); yn(nm,2)=y(j,2); zn(nm,2)=z(j,2); efln(nm,2)=efl(j,2);
   xn(nm,3)=x23;    yn(nm,3)=y23;    zn(nm,3)=z23;    efln(nm,3)=efl23;
   xn(nm,4)=x24;    yn(nm,4)=y24;    zn(nm,4)=z24;    efln(nm,4)=efl24;

   nm = nm+1;  % third subelement

   xn(nm,1)=x14;    yn(nm,1)=y14;    zn(nm,1)=z14;    efln(nm,1)=efl14;
   xn(nm,2)=x24;    yn(nm,2)=y24;    zn(nm,2)=z24;    efln(nm,2)=efl24;
   xn(nm,3)=x34;    yn(nm,3)=y34;    zn(nm,3)=z34;    efln(nm,3)=efl34;
   xn(nm,4)=x(j,4); yn(nm,4)=y(j,4); zn(nm,4)=z(j,4); efln(nm,4)=efl(j,4);

   nm = nm+1;  % fourth subelement

   xn(nm,1)=x(j,1); yn(nm,1)=y(j,1); zn(nm,1)=z(j,1); efln(nm,1)=efl(j,1);
   xn(nm,2)=x12;    yn(nm,2)=y12;    zn(nm,2)=z12;    efln(nm,2)=efl12;
   xn(nm,3)=x13;    yn(nm,3)=y13;    zn(nm,3)=z13;    efln(nm,3)=efl13;
   xn(nm,4)=x14;    yn(nm,4)=y14;    zn(nm,4)=z14;    efln(nm,4)=efl14;

   nm = nm+1;  % fifth subelement

   xn(nm,1)=x34; yn(nm,1)=y34; zn(nm,1)=z34; efln(nm,1)=efl34;
   xn(nm,2)=x24; yn(nm,2)=y24; zn(nm,2)=z24; efln(nm,2)=efl24;
   xn(nm,3)=x14; yn(nm,3)=y14; zn(nm,3)=z14; efln(nm,3)=efl14;
   xn(nm,4)=x13; yn(nm,4)=y13; zn(nm,4)=z13; efln(nm,4)=efl13;

   nm = nm+1;  % sixth subelement

   xn(nm,1)=x14; yn(nm,1)=y14; zn(nm,1)=z14; efln(nm,1)=efl14;
   xn(nm,2)=x12; yn(nm,2)=y12; zn(nm,2)=z12; efln(nm,2)=efl12;
   xn(nm,3)=x13; yn(nm,3)=y13; zn(nm,3)=z13; efln(nm,3)=efl13;
   xn(nm,4)=x24; yn(nm,4)=y24; zn(nm,4)=z24; efln(nm,4)=efl24;

   nm = nm+1;  % seventh subelement

   xn(nm,1)=x12; yn(nm,1)=y12; zn(nm,1)=z12; efln(nm,1)=efl12;
   xn(nm,2)=x23; yn(nm,2)=y23; zn(nm,2)=z23; efln(nm,2)=efl23;
   xn(nm,3)=x13; yn(nm,3)=y13; zn(nm,3)=z13; efln(nm,3)=efl13;
   xn(nm,4)=x24; yn(nm,4)=y24; zn(nm,4)=z24; efln(nm,4)=efl24;

   nm = nm+1;  % eighth subelement

   xn(nm,1)=x13; yn(nm,1)=y13; zn(nm,1)=z13; efln(nm,1)=efl13;
   xn(nm,2)=x24; yn(nm,2)=y24; zn(nm,2)=z24; efln(nm,2)=efl24;
   xn(nm,3)=x23; yn(nm,3)=y23; zn(nm,3)=z23; efln(nm,3)=efl23;
   xn(nm,4)=x34; yn(nm,4)=y34; zn(nm,4)=z34; efln(nm,4)=efl34;

 end % end of loop over current elements

%---

 ne = 8*ne;  % number of elements has increased
             % by a factor of eight

 for k=1:ne     % relabel the new points
                % and put them in the master list
                % project boundary nodes onto the unit sphere
    for l=1:4
       x(k,l) =   xn(k,l);
       y(k,l) =   yn(k,l);
       z(k,l) =   zn(k,l);
     efl(k,l) = efln(k,l);
      if(efl(k,l) == 1)
        rad = sqrt(x(k,l)^2+y(k,l)^2+z(k,l)^2);
        fc = act/rad;
        x(k,l) = fc*x(k,l);
        y(k,l) = fc*y(k,l);
        z(k,l) = fc*z(k,l);
      end
    end

 end

end % end do of refinement loop
end % end if of refinement loop

%---------------------------------------------------
% define the global nodes and the connectivity table
%---------------------------------------------------

% four nodes of the first element are entered mannualy

p(1,1)=x(1,1); p(1,2)=y(1,1); p(1,3)=z(1,1); gfl(1)=efl(1,1);
p(2,1)=x(1,2); p(2,2)=y(1,2); p(2,3)=z(1,2); gfl(2)=efl(1,2);
p(3,1)=x(1,3); p(3,2)=y(1,3); p(3,3)=z(1,3); gfl(3)=efl(1,3);
p(4,1)=x(1,4); p(4,2)=y(1,4); p(4,3)=z(1,4); gfl(4)=efl(1,4);

c(1,1) = 1;  % first  node of first element is global node 1
c(1,2) = 2;  % second node of first element is global node 2
c(1,3) = 3;  % third  node of first element is global node 3
c(1,4) = 4;  % fourth node of first element is global node 4

ng = 4;

%---
% loop over further elements
% Iflag=0 will signal a new global node
%---

eps = 0.00001;  % tolerance

for i=2:ne        % loop over the rest of the elements
 for j=1:4          % loop over element nodes

 Iflag=0;

 for k=1:ng

  if(abs(x(i,j)-p(k,1)) < eps)
   if(abs(y(i,j)-p(k,2)) < eps)
    if(abs(z(i,j)-p(k,3)) < eps)

     Iflag = 1;    % the node has been recorded previously
     c(i,j) = k;   % the jth local node of element i is the kth global node

    end
   end
  end

 end

 if(Iflag==0)  % record the node
   ng = ng+1;
   p(ng,1) =   x(i,j);
   p(ng,2) =   y(i,j);
   p(ng,3) =   z(i,j);
   gfl(ng) = efl(i,j);
   c(i,j) = ng;   % the jth local node of element is the new global node
 end

 end
end  % end of loop over elements

%-----
% done
%-----

return;

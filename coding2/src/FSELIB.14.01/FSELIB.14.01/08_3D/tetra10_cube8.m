function [ne,ng,p,c,efl,gfl] = tetra10_cube8 (ndiv)

%=====================================================
% Discretize a cube into 10-node tetrahedral elements
% by successive subdivision of an 12-element structure
% (root set) into 8 descendant elements
% 
% LEGEND:
%
% efl: element flag
% gfl: global node flag
%====================================

%---------------------------------------
% parental structure with eight elements
%---------------------------------------

ne=12;

x(1,1)= 0.0; y(1,1)= 0.0; z(1,1)= 0.0; efl(1,1)=0;  % first element
x(1,2)= 1.0; y(1,2)=-1.0; z(1,2)= 1.0; efl(1,2)=1;
x(1,3)= 1.0; y(1,3)= 1.0; z(1,3)=-1.0; efl(1,3)=1;
x(1,4)= 1.0; y(1,4)= 1.0; z(1,4)= 1.0; efl(1,4)=1;

x(2,1)= 0.0; y(2,1)= 0.0; z(2,1)= 0.0; efl(2,1)=0;  % second element
x(2,2)= 1.0; y(2,2)=-1.0; z(2,2)= 1.0; efl(2,2)=1;
x(2,3)= 1.0; y(2,3)=-1.0; z(2,3)=-1.0; efl(2,3)=1;
x(2,4)= 1.0; y(2,4)= 1.0; z(2,4)=-1.0; efl(2,4)=1;

x(3,1)= 0.0; y(3,1)= 0.0; z(3,1)= 0.0; efl(3,1)=0;  % third element
x(3,2)= 1.0; y(3,2)=-1.0; z(3,2)= 1.0; efl(3,2)=1;
x(3,3)=-1.0; y(3,3)=-1.0; z(3,3)= 1.0; efl(3,3)=1;
x(3,4)= 1.0; y(3,4)=-1.0; z(3,4)=-1.0; efl(3,4)=1;

x(4,1)= 0.0; y(4,1)= 0.0; z(4,1)= 0.0; efl(4,1)=0;  % fourth element
x(4,2)=-1.0; y(4,2)=-1.0; z(4,2)= 1.0; efl(4,2)=1;
x(4,3)=-1.0; y(4,3)=-1.0; z(4,3)=-1.0; efl(4,3)=1;
x(4,4)= 1.0; y(4,4)=-1.0; z(4,4)=-1.0; efl(4,4)=1;

x(5,1)= 0.0; y(5,1)= 0.0; z(5,1)= 0.0; efl(5,1)=0;  % fifth element
x(5,2)=-1.0; y(5,2)= 1.0; z(5,2)= 1.0; efl(5,2)=1;
x(5,3)= 1.0; y(5,3)= 1.0; z(5,3)= 1.0; efl(5,3)=1;
x(5,4)= 1.0; y(5,4)= 1.0; z(5,4)=-1.0; efl(5,4)=1;

x(6,1)= 0.0; y(6,1)= 0.0; z(6,1)= 0.0; efl(6,1)=0;  % sixth element
x(6,2)=-1.0; y(6,2)= 1.0; z(6,2)=-1.0; efl(6,2)=1;
x(6,3)=-1.0; y(6,3)= 1.0; z(6,3)= 1.0; efl(6,3)=1;
x(6,4)= 1.0; y(6,4)= 1.0; z(6,4)=-1.0; efl(6,4)=1;

x(7,1)= 0.0; y(7,1)= 0.0; z(7,1)= 0.0; efl(7,1)=0;  % seventh element
x(7,2)=-1.0; y(7,2)= 1.0; z(7,2)=-1.0; efl(7,2)=1;
x(7,3)=-1.0; y(7,3)=-1.0; z(7,3)= 1.0; efl(7,3)=1;
x(7,4)=-1.0; y(7,4)= 1.0; z(7,4)= 1.0; efl(7,4)=1;

x(8,1)= 0.0; y(8,1)= 0.0; z(8,1)= 0.0; efl(8,1)=0;  % eighth element
x(8,2)=-1.0; y(8,2)=-1.0; z(8,2)=-1.0; efl(8,2)=1;
x(8,3)=-1.0; y(8,3)=-1.0; z(8,3)= 1.0; efl(8,3)=1;
x(8,4)=-1.0; y(8,4)= 1.0; z(8,4)=-1.0; efl(8,4)=1;

x(9,1)= 0.0; y(9,1)= 0.0; z(9,1)= 0.0; efl(9,1)=0;  % ninth element
x(9,2)=-1.0; y(9,2)=-1.0; z(9,2)= 1.0; efl(9,2)=1;
x(9,3)= 1.0; y(9,3)=-1.0; z(9,3)= 1.0; efl(9,3)=1;
x(9,4)=-1.0; y(9,4)= 1.0; z(9,4)= 1.0; efl(9,4)=1;

x(10,1)= 0.0; y(10,1)= 0.0; z(10,1)= 0.0; efl(10,1)=0;  % tenth element
x(10,2)=-1.0; y(10,2)= 1.0; z(10,2)= 1.0; efl(10,2)=1;
x(10,3)= 1.0; y(10,3)=-1.0; z(10,3)= 1.0; efl(10,3)=1;
x(10,4)= 1.0; y(10,4)= 1.0; z(10,4)= 1.0; efl(10,4)=1;

x(11,1)= 0.0; y(11,1)= 0.0; z(11,1)= 0.0; efl(11,1)=0;  % eleventh element
x(11,2)= 1.0; y(11,2)=-1.0; z(11,2)=-1.0; efl(11,2)=1;
x(11,3)=-1.0; y(11,3)=-1.0; z(11,3)=-1.0; efl(11,3)=1;
x(11,4)=-1.0; y(11,4)= 1.0; z(11,4)=-1.0; efl(11,4)=1;

x(12,1)= 0.0; y(12,1)= 0.0; z(12,1)= 0.0; efl(12,1)=0;  % twelfth element
x(12,2)= 1.0; y(12,2)=-1.0; z(12,2)=-1.0; efl(12,2)=1;
x(12,3)=-1.0; y(12,3)= 1.0; z(12,3)=-1.0; efl(12,3)=1;
x(12,4)= 1.0; y(12,4)= 1.0; z(12,4)=-1.0; efl(12,4)=1;

%------------------
% element midpoints
%------------------

for i=1:ne

 x(i,5)  = 0.5*(x(i,1)+x(i,2));
 y(i,5)  = 0.5*(y(i,1)+y(i,2));
 z(i,5)  = 0.5*(z(i,1)+z(i,2));
 efl(i,5)  = 0;

 x(i,6)  = 0.5*(x(i,2)+x(i,3));
 y(i,6)  = 0.5*(y(i,2)+y(i,3));
 z(i,6)  = 0.5*(z(i,2)+z(i,3));
 efl(i,6)  = 1;

 x(i,7)  = 0.5*(x(i,1)+x(i,3));
 y(i,7)  = 0.5*(y(i,1)+y(i,3));
 z(i,7)  = 0.5*(z(i,1)+z(i,3));
 efl(i,7)  = 0;

 x(i,8)  = 0.5*(x(i,1)+x(i,4));
 y(i,8)  = 0.5*(y(i,1)+y(i,4));
 z(i,8)  = 0.5*(z(i,1)+z(i,4));
 efl(i,8)  = 0;

 x(i,9)  = 0.5*(x(i,2)+x(i,4));
 y(i,9)  = 0.5*(y(i,2)+y(i,4));
 z(i,9)  = 0.5*(z(i,2)+z(i,4));
 efl(i,9)  = 1;

 x(i,10) = 0.5*(x(i,3)+x(i,4));
 y(i,10) = 0.5*(y(i,3)+y(i,4));
 z(i,10) = 0.5*(z(i,3)+z(i,4));
 efl(i,10) = 1;

end

%----------------
% refinement loop
%----------------

if(ndiv > 0)

for i=1:ndiv

 nm = 0; % count the new elements arising by each refinement loop
         % eight elements will be generated in each pass

 for j=1:ne   % loop over current elements

  % assign vertex nodes to sub-elements
  % these will become the "new" elements

   nm = nm+1;  % first sub-element

   xn(nm,1)=x(j,7);  yn(nm,1)=y(j,7);  zn(nm,1)=z(j,7);
                                   efln(nm,1)=efl(j,7);
   xn(nm,2)=x(j,6);  yn(nm,2)=y(j,6);  zn(nm,2)=z(j,6);
                                   efln(nm,2)=efl(j,6);
   xn(nm,3)=x(j,3);  yn(nm,3)=y(j,3);  zn(nm,3)=z(j,3);
                                   efln(nm,3)=efl(j,3);
   xn(nm,4)=x(j,10); yn(nm,4)=y(j,10); zn(nm,4)=z(j,10);
                                   efln(nm,4)=efl(j,10);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % second subelement

   xn(nm,1)=x(j,5);  yn(nm,1)=y(j,5); zn(nm,1)=z(j,5);
                                  efln(nm,1)=efl(j,5);
   xn(nm,2)=x(j,2);  yn(nm,2)=y(j,2); zn(nm,2)=z(j,2);
                                  efln(nm,2)=efl(j,2);
   xn(nm,3)=x(j,6);  yn(nm,3)=y(j,6); zn(nm,3)=z(j,6);
                                  efln(nm,3)=efl(j,6);
   xn(nm,4)=x(j,9);  yn(nm,4)=y(j,9); zn(nm,4)=z(j,9);
                                  efln(nm,4)=efl(j,9);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % third subelement

   xn(nm,1)=x(j,8);  yn(nm,1)=y(j,8);  zn(nm,1)=z(j,8);
                                   efln(nm,1)=efl(j,8);
   xn(nm,2)=x(j,9);  yn(nm,2)=y(j,9);  zn(nm,2)=z(j,9);
                                   efln(nm,2)=efl(j,9);
   xn(nm,3)=x(j,10); yn(nm,3)=y(j,10); zn(nm,3)=z(j,10);
                                   efln(nm,3)=efl(j,10);
   xn(nm,4)=x(j,4);  yn(nm,4)=y(j,4);  zn(nm,4)=z(j,4);
                                   efln(nm,4)=efl(j,4);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % fourth subelement

   xn(nm,1)=x(j,1); yn(nm,1)=y(j,1); zn(nm,1)=z(j,1);
                                 efln(nm,1)=efl(j,1);
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); zn(nm,2)=z(j,5);
                                 efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,7); yn(nm,3)=y(j,7); zn(nm,3)=z(j,7);
                                 efln(nm,3)=efl(j,7);
   xn(nm,4)=x(j,8); yn(nm,4)=y(j,8); zn(nm,4)=z(j,8);
                                 efln(nm,4)=efl(j,8);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % fifth subelement

   xn(nm,1)=x(j,10); yn(nm,1)=y(j,10); zn(nm,1)=z(j,10);
                                   efln(nm,1)=efl(j,10);
   xn(nm,2)=x(j,9);  yn(nm,2)=y(j,9);  zn(nm,2)=z(j,9);
                                   efln(nm,2)=efl(j,9);
   xn(nm,3)=x(j,8);  yn(nm,3)=y(j,8);  zn(nm,3)=z(j,8);
                                   efln(nm,3)=efl(j,8);
   xn(nm,4)=x(j,7);  yn(nm,4)=y(j,7);  zn(nm,4)=z(j,7);
                                   efln(nm,4)=efl(j,7);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % sixth subelement

   xn(nm,1)=x(j,8); yn(nm,1)=y(j,8); zn(nm,1)=z(j,8);
                                 efln(nm,1)=efl(j,8);
   xn(nm,2)=x(j,5); yn(nm,2)=y(j,5); zn(nm,2)=z(j,5);
                                 efln(nm,2)=efl(j,5);
   xn(nm,3)=x(j,7); yn(nm,3)=y(j,7); zn(nm,3)=z(j,7);
                                 efln(nm,3)=efl(j,7);
   xn(nm,4)=x(j,9); yn(nm,4)=y(j,9); zn(nm,4)=z(j,9);
                                 efln(nm,4)=efl(j,9);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % seventh subelement

   xn(nm,1)=x(j,5); yn(nm,1)=y(j,5); zn(nm,1)=z(j,5);
                                 efln(nm,1)=efl(j,5);
   xn(nm,2)=x(j,6); yn(nm,2)=y(j,6); zn(nm,2)=z(j,6);
                                 efln(nm,2)=efl(j,6);
   xn(nm,3)=x(j,7); yn(nm,3)=y(j,7); zn(nm,3)=z(j,7);
                                 efln(nm,3)=efl(j,7);
   xn(nm,4)=x(j,9); yn(nm,4)=y(j,9); zn(nm,4)=z(j,9);
                                 efln(nm,4)=efl(j,9);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

   nm = nm+1;  % eighth subelement

   xn(nm,1)=x(j,7);  yn(nm,1)=y(j,7);  zn(nm,1)=z(j,7);
                                   efln(nm,1)=efl(j,7);
   xn(nm,2)=x(j,9);  yn(nm,2)=y(j,9);  zn(nm,2)=z(j,9);
                                   efln(nm,2)=efl(j,9);
   xn(nm,3)=x(j,6);  yn(nm,3)=y(j,6);  zn(nm,3)=z(j,6);
                                   efln(nm,3)=efl(j,6);
   xn(nm,4)=x(j,10); yn(nm,4)=y(j,10); zn(nm,4)=z(j,10);
                                   efln(nm,4)=efl(j,10);

   xn(nm,5)  = 0.5*(xn(nm,1)+xn(nm,2));
   yn(nm,5)  = 0.5*(yn(nm,1)+yn(nm,2));
   zn(nm,5)  = 0.5*(zn(nm,1)+zn(nm,2)); efln(nm,5) = 0;

   xn(nm,6)  = 0.5*(xn(nm,2)+xn(nm,3));
   yn(nm,6)  = 0.5*(yn(nm,2)+yn(nm,3));
   zn(nm,6)  = 0.5*(zn(nm,2)+zn(nm,3)); enfln(nm,6) = 0;

   xn(nm,7)  = 0.5*(xn(nm,1)+xn(nm,3));
   yn(nm,7)  = 0.5*(yn(nm,1)+yn(nm,3));
   zn(nm,7)  = 0.5*(zn(nm,1)+zn(nm,3)); efln(nm,7) = 0;

   xn(nm,8)  = 0.5*(xn(nm,1)+xn(nm,4));
   yn(nm,8)  = 0.5*(yn(nm,1)+yn(nm,4));
   zn(nm,8)  = 0.5*(zn(nm,1)+zn(nm,4)); efln(nm,8) = 0;

   xn(nm,9)  = 0.5*(xn(nm,2)+xn(nm,4));
   yn(nm,9)  = 0.5*(yn(nm,2)+yn(nm,4));
   zn(nm,9)  = 0.5*(zn(nm,2)+zn(nm,4)); efln(nm,9) = 0;

   xn(nm,10) = 0.5*(xn(nm,3)+xn(nm,4));
   yn(nm,10) = 0.5*(yn(nm,3)+yn(nm,4));
   zn(nm,10) = 0.5*(zn(nm,3)+zn(nm,4)); efln(nm,10) = 0;

   if(efln(nm,1)==1 & efln(nm,2)==1) efln(nm,5)  = 1; end
   if(efln(nm,2)==1 & efln(nm,3)==1) efln(nm,6)  = 1; end
   if(efln(nm,1)==1 & efln(nm,3)==1) efln(nm,7)  = 1; end
   if(efln(nm,1)==1 & efln(nm,4)==1) efln(nm,8)  = 1; end
   if(efln(nm,2)==1 & efln(nm,4)==1) efln(nm,9)  = 1; end
   if(efln(nm,3)==1 & efln(nm,4)==1) efln(nm,10) = 1; end

 end % end of loop over current elements

%---

 ne = 8*ne;  % number of elements has increased
             % by a factor of eight

 for k=1:ne     % relabel the new points
                % and put them in the master list
                % project boundary nodes onto a sphere
                % of radius "act"
    for l=1:10
       x(k,l) =   xn(k,l);
       y(k,l) =   yn(k,l);
       z(k,l) =   zn(k,l);
     efl(k,l) = efln(k,l);
    end

  end

end % end do of refinement loop
end % end if of refinement loop

%---------------------------------------------------
% define the global nodes and the connectivity table
%---------------------------------------------------

% ten nodes of the first element are entered mannualy

p(1,1) =x(1,1);  p(1,2) =y(1,1);  p(1,3) =z(1,1);  gfl(1)=efl(1,1);
p(2,1) =x(1,2);  p(2,2) =y(1,2);  p(2,3) =z(1,2);  gfl(2)=efl(1,2);
p(3,1) =x(1,3);  p(3,2) =y(1,3);  p(3,3) =z(1,3);  gfl(3)=efl(1,3);
p(4,1) =x(1,4);  p(4,2) =y(1,4);  p(4,3) =z(1,4);  gfl(4)=efl(1,4);
p(5,1) =x(1,5);  p(5,2) =y(1,5);  p(5,3) =z(1,5);  gfl(5)=efl(1,5);
p(6,1) =x(1,6);  p(6,2) =y(1,6);  p(6,3) =z(1,6);  gfl(6)=efl(1,6);
p(7,1) =x(1,7);  p(7,2) =y(1,7);  p(7,3) =z(1,7);  gfl(7)=efl(1,7);
p(8,1) =x(1,8);  p(8,2) =y(1,8);  p(8,3) =z(1,8);  gfl(8)=efl(1,8);
p(9,1) =x(1,9);  p(9,2) =y(1,9);  p(9,3) =z(1,9);  gfl(9)=efl(1,9);
p(10,1)=x(1,10); p(10,2)=y(1,10); p(10,3)=z(1,10); gfl(10)=efl(1,10);

c(1,1)  = 1;  % first  node of first element is global node 1
c(1,2)  = 2;  % second node of first element is global node 2
c(1,3)  = 3;  % third  node of first element is global node 3
c(1,4)  = 4;  % fourth node of first element is global node 4
c(1,5)  = 5;  %
c(1,6)  = 6;  %
c(1,7)  = 7;  %
c(1,8)  = 8;  % 
c(1,9)  = 9;  % 
c(1,10) = 10; % 

ng = 10;

%---
% loop over further elements
% Iflag=0 will signal a new global node
%---

eps = 0.00001;  % tolerance

for i=2:ne        % loop over the rest of the elements
 for j=1:10          % loop over element nodes

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
   gfl(ng)   = efl(i,j);
   c(i,j) = ng;   % the jth local node of element is the new global node
 end

 end
end  % end of loop over elements

%-----
% done
%-----

return;

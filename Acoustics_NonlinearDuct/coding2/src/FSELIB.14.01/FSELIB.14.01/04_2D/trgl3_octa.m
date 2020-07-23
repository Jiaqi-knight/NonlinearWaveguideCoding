function [Npts,Nelm,p,ne,n,nbe] = trgl3_octa (Ndiv)

%=========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%=========================================

%--------------------------------------------
% Triangulation of the unit sphere
% into three -node triangles
% by the successive subdivision of a regular octahedron
%
% SYMBOLS:
% -------
%
%  Ndiv .... level of discretization of an octahedron
%            nvid = 0 gives 8 elements
%
%  Npts .... number of nodes
%  Nelm .... number of surface elements
%
%  x(i,j), y(i,j), z(i,j) .... 
%
%    Cartesian coordinates of local node j
%    on element i, j = 1,2,3 and i = 1,...,Nelm
%
%  p(i,j) .... 
%
%  Cartesian coordinates of global node i
%  where j=1,2,3, with: x = p(i,1), y = p(i,2), z = p(i,3)
%                                 
%
%  n(i,j) .... global node number of local node number j on element i,
%              where j=1,2,3
%
%  ne(i,j) ... ne(i,1) is the number of elements touching global node i.
%              ne(i,2:ne(i,1)) are the corresponding element labels 
%
%  nbe(i,j) .. label of element sharing side j of element i
%              where j = 1, 2, 3
%--------------------------------------------------------------------

%----------------------------------------
% Begin with the zeroth-level
% discretization (8 elements)
%
% Nodes are set manually on the unit sphere
%----------------------------------------

      Nelm = 8;

%---
%  vertex nodes in the upper half of xz plane
%---

      x(1,1)= 0.0D0;   % first element
      y(1,1)= 0.0D0;
      z(1,1)= 1.0D0;

      x(1,2)= 1.0D0; 
      y(1,2)= 0.0D0;
      z(1,2)= 0.0D0;

      x(1,3)= 0.0D0;
      y(1,3)= 1.0D0;
      z(1,3)= 0.0D0;
%---
      x(5,1)= 1.0D0;   %  fifth element
      y(5,1)= 0.0D0;
      z(5,1)= 0.0D0;

      x(5,2)= 0.0D0;
      y(5,2)= 0.0D0;
      z(5,2)=-1.0D0;

      x(5,3)= 0.0D0;
      y(5,3)= 1.0D0;
      z(5,3)= 0.0D0;
%---
      x(6,1)= 0.0D0;  % sixth element
      y(6,1)= 0.0D0;
      z(6,1)=-1.0D0;

      x(6,2)=-1.0D0;
      y(6,2)= 0.0D0;
      z(6,2)= 0.0D0;

      x(6,3)= 0.0D0;
      y(6,3)= 1.0D0;
      z(6,3)= 0.0D0;
%---
      x(2,1)=-1.0D0;    % second element
      y(2,1)= 0.0D0;
      z(2,1)= 0.0D0;

      x(2,2)= 0.0D0;
      y(2,2)= 0.0D0;
      z(2,2)= 1.0D0;

      x(2,3)= 0.0D0;
      y(2,3)= 1.0D0;
      z(2,3)= 0.0D0;

%---
%  vertex nodes in the lower half xz plane
%---

      x(4,1)= 0.0D0;    % fourth element
      y(4,1)= 0.0D0;
      z(4,1)= 1.0D0;

      x(4,2)= 0.0D0;
      y(4,2)=-1.0D0;
      z(4,2)= 0.0D0;

      x(4,3)= 1.0D0;
      y(4,3)= 0.0D0;
      z(4,3)= 0.0D0;

%---

      x(8,1)= 1.0D0 ;   % eighth element
      y(8,1)= 0.0D0;
      z(8,1)= 0.0D0;

      x(8,2)= 0.0D0;
      y(8,2)=-1.0D0;
      z(8,2)= 0.0D0;

      x(8,3)= 0.0D0;
      y(8,3)= 0.0D0;
      z(8,3)=-1.0D0;

%---

      x(7,1)= 0.0D0;   % seventh element
      y(7,1)= 0.0D0;
      z(7,1)=-1.0D0;

      x(7,2)= 0.0D0;
      y(7,2)=-1.0D0;
      z(7,2)= 0.0D0;

      x(7,3)=-1.0D0;
      y(7,3)= 0.0D0;
      z(7,3)= 0.0D0;

%---

      x(3,1)=-1.0D0;    % third element
      y(3,1)= 0.0D0;
      z(3,1)= 0.0D0;

      x(3,2)= 0.0D0;
      y(3,2)=-1.0D0;
      z(3,2)= 0.0D0;

      x(3,3)= 0.0D0;
      y(3,3)= 0.0D0;
      z(3,3)= 1.0D0;

%------------------------------------------
% Compute the mid-points of the three sides
% of the 8 first-generation elements
%
% mid-points are numbered 4, 5, 6
%------------------------------------------

      for i=1:Nelm

       x(i,4) = 0.5D0*(x(i,1)+x(i,2));
       y(i,4) = 0.5D0*(y(i,1)+y(i,2));
       z(i,4) = 0.5D0*(z(i,1)+z(i,2));

       x(i,5) = 0.5D0*(x(i,2)+x(i,3));
       y(i,5) = 0.5D0*(y(i,2)+y(i,3));
       z(i,5) = 0.5D0*(z(i,2)+z(i,3));

       x(i,6) = 0.5D0*(x(i,3)+x(i,1));
       y(i,6) = 0.5D0*(y(i,3)+y(i,1));
       z(i,6) = 0.5D0*(z(i,3)+z(i,1));

      end

%---
% project the nodes onto the unit sphere
%---

      for k=1:Nelm
        for l=1:6
         rad = sqrt(x(k,l)^2+y(k,l)^2+z(k,l)^2); 
         x(k,l) = x(k,l)/rad;
         y(k,l) = y(k,l)/rad;
         z(k,l) = z(k,l)/rad;
       end
      end

%-------------------------------------------
if(Ndiv>0)
%-------------------------------------------

%-------------------------------------------
% Compute the local element node coordinates
% for discretization levels 1 through Ndiv
%-------------------------------------------

      for i=1:Ndiv     % loop over refinement levels

       nm = 0;        % counts the new elements arising in each refinement loop
                      % four elements will be generated in each pass
       for j=1:Nelm    % loop over old elements

%---
% assign corner points to sub-elements
% these will become the "new" elements
%---

        nm = nm+1;

        xn(nm,1) = x(j,1);     %  first sub-element
        yn(nm,1) = y(j,1);
        zn(nm,1) = z(j,1);

        xn(nm,2) = x(j,4);
        yn(nm,2) = y(j,4);
        zn(nm,2) = z(j,4);

        xn(nm,3) = x(j,6);
        yn(nm,3) = y(j,6);
        zn(nm,3) = z(j,6);

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2));
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2));
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2));

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3));
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3));
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3));

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1));
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1));
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1));

        nm = nm+1;

        xn(nm,1) = x(j,4);      %  second sub-element
        yn(nm,1) = y(j,4);
        zn(nm,1) = z(j,4);

        xn(nm,2) = x(j,2);
        yn(nm,2) = y(j,2);
        zn(nm,2) = z(j,2);

        xn(nm,3) = x(j,5);
        yn(nm,3) = y(j,5);
        zn(nm,3) = z(j,5);

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2));
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2));
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2));

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3));
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3));
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3));

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1));
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1));
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1));

        nm = nm+1;

        xn(nm,1) = x(j,6);     %  third sub-element
        yn(nm,1) = y(j,6);
        zn(nm,1) = z(j,6);
 
        xn(nm,2) = x(j,5);
        yn(nm,2) = y(j,5);
        zn(nm,2) = z(j,5);

        xn(nm,3) = x(j,3);
        yn(nm,3) = y(j,3);
        zn(nm,3) = z(j,3);

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2));
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2));
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2));

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3));
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3));
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3));

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1));
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1));
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1));

        nm = nm+1;

        xn(nm,1) = x(j,4);     %  fourth sub-element
        yn(nm,1) = y(j,4);
        zn(nm,1) = z(j,4);

        xn(nm,2) = x(j,5);
        yn(nm,2) = y(j,5);
        zn(nm,2) = z(j,5);

        xn(nm,3) = x(j,6);
        yn(nm,3) = y(j,6);
        zn(nm,3) = z(j,6);

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2));    % mid points
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2));
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2));

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3));
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3));
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3));

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1));
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1));
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1));

       end %  of loop over old elements

%--------------------------------------
% number of elements has been increased
% by a factor of 4
%--------------------------------------

       Nelm = 4*Nelm;

%-----------------------------------
% relabel the new points (xn -> x)
% and place them in the master list
%----------------------------------

       for k=1:Nelm
        for l=1:6

         x(k,l) = xn(k,l);
         y(k,l) = yn(k,l);
         z(k,l) = zn(k,l);

%--- project onto the sphere

         rad = sqrt(x(k,l)^2+y(k,l)^2+z(k,l)^2);
         x(k,l) = x(k,l)/rad;
         y(k,l) = y(k,l)/rad;
         z(k,l) = z(k,l)/rad;

         xn(k,l) = 0.0D0;   % zero just in case
         yn(k,l) = 0.0D0;
         zn(k,l) = 0.0D0;

        end
       end

%----------

      end %  of refinement loop

%-----------------------------------------
end
%-----------------------------------------

%-----------------------------------
% Generate a list of global nodes by looping 
% over all elements
% and adding nodes not found in the list.
%
% Fill in the connectivity table n(i,j) 
% containing node numbers of element points 1,2,3
%-----------------------------------

%---
% three nodes of the first element are
% entered mannualy
%---

      p(1,1) = x(1,1);
      p(1,2) = y(1,1);
      p(1,3) = z(1,1);

      p(2,1) = x(1,2);
      p(2,2) = y(1,2);
      p(2,3) = z(1,2);

      p(3,1) = x(1,3);
      p(3,2) = y(1,3);
      p(3,3) = z(1,3);

      n(1,1) = 1;  % first  node of first element is global node 1
      n(1,2) = 2;  % second node of first element is global node 2
      n(1,3) = 3;  % third  node of first element is global node 3

      Npts = 3;

%---
% loop over further elements
%
% Iflag=0 will signal a new global node
%
% n(i,j): global label on jth node on ith element
%---

      for i=2:Nelm     % loop over elements
       for j=1:3        % loop over element nodes

        Iflag=0;

         for k=1:Npts
          if(abs(x(i,j)-p(k,1)) <= eps)
           if(abs(y(i,j)-p(k,2)) <= eps) 
            if(abs(z(i,j)-p(k,3)) <= eps)

             Iflag = 1;    % the node has been recorded previously
             n(i,j) = k;   % the jth local node of element i
                          % is the kth global node 
            end
           end
          end
         end
        
         if(Iflag == 0)  % record the node

          Npts = Npts+1;         % one more global node

          p(Npts,1) = x(i,j);
          p(Npts,2) = y(i,j);
          p(Npts,3) = z(i,j);

          n(i,j) = Npts;   % the jth local node of element i
                          % is the new global node 
         end

       end
      end    %  of loop over elements

%----------------------------------
% Generate connectivity table: ne(i,j)
% for elements touching global node i
%
% ne(i,1) is the number of elements touching  global node i
% ne(i,j) for j=2, ..., ne(i,1)+1
%         are the corresponding element labels
%----------------------------------

%---
% initialize
%---

      for i=1:Npts
       for j=1:3
        ne(i,j) = 0;
       end
      end

%---
% loop over global nodes
%---

      for i=1:Npts 

       ne(i,1) = 0;
       Icount = 1;

       for j=1:Nelm     % loop over elements
        for k=1:3       % loop over element nodes

           if(n(j,k)==i)
             ne(i,1) = ne(i,1)+1;
             Icount =Icount+1;
             ne(i,Icount) = j;
           end

        end
       end
     
      end  %  of loop over global nodes

%------------------------------------------
% Generate connectivity table nbe(i,j) for 
%
%  nbe(i,j) .. label of element sharing side j of element i
%              where j = 1, 2, 3
%------------------------------------------

%---
% initialize
%---

      for i=1:Nelm
       for j=1:3
        nbe(i,j) = 0;
       end
      end

%---
% loop over elements
%---

      for i=1:Nelm            %  loop over elements
       jcount = 1;
       for j=4:6              %  loop over mid-points

        for k=1:Nelm           %  test element
         if(k~=i) % not a self-element
          for l=4:6             %  loop over mid-points

          if(abs(x(i,j)-x(k,l)) <= eps) 
           if(abs(y(i,j)-y(k,l)) <= eps) 
            if(abs(z(i,j)-z(k,l)) <= eps)
             nbe(i,jcount) = k;
            end
           end
          end

          end
         end
        end            %  end of test element

        if(nbe(i,jcount)~=0) 
         jcount = jcount+1;
        end

       end
      end  %  of loop over elements

%-------------------------------------------
% project points p(i,j) onto the unit sphere
%-------------------------------------------

      for i=1:Npts

       rad = sqrt(p(i,1)^2+p(i,2)^2+p(i,3)^2);

       p(i,1) = p(i,1)/rad;
       p(i,2) = p(i,2)/rad;
       p(i,3) = p(i,3)/rad;

      end
%-----
% Done
%-----

return

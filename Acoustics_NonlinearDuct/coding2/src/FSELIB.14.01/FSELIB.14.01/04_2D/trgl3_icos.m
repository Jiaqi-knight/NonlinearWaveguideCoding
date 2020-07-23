function [Npts,Nelm,p,ne,n,nbe] = trgl3_icos (Ndiv)

%-----------------------------------------
% FDLIB BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%----------------------------------------

%--------------------------------------------
% Triangulation of the unit sphere
% into three-node quadratic triangles
% by subdividing a regular icosahedron
%
% SYMBOLS:
% -------
%
%  Ndiv .... level of discretization of icosahedron
%            Nvid = 0 gives 20 elements
%
%  Npts .... number of nodes
%  Nelm .... number of surface elements
%
%  x(i,j), y(i,j), z(i,j) .... 
%
%        Cartesian coordinates of local node j (1,2,3)
%        on element i (1...Nelm)
%
%  p(i,j) .... Cartesian coordinates of global node i
%              where j=1,2,3, with
%              x = p(i,1)  y = p(i,2)  z = p(i,3)
%                           
%  n(i,j) .... global node label of local node number j on element i,
%              where j=1,2,3
%
%  ne(i,j) ... ne(i,1) is the number of elements touching global node i.
%              ne(i, 2:ne(i,1)+1) are the corresponding element labels 
%
%  nbe(i,j) .. label of element sharing side j of element i
%              where j = 1, 2, 3
%-------------------------------------------------------------

%----------------------------------------
% Begin with the zeroth-level
% discretization (20 elements)
%
% Nodes are set manually on the unit sphere
%----------------------------------------

%---
% the icosahedron has 12 vertices
%---

  ru = 0.25D0*sqrt(10.0D0+2.0D0*sqrt(5.0D0));
  rm = 0.25D0*(1.0D0+sqrt(5.0D0));   % twice the golden ratio
  c0 = 0.0;
  c1 = ru;
  c2 = 2.0*ru/sqrt(5.0D0);
  c3 = rm;
  c4 = ru/sqrt(5.0D0);
  c5 = sqrt( ru^2-c3^2-c4^2);
  c6 = 0.5D0;
  c7 = sqrt( ru^2-c4^2-c6^2);

      VX(1) =  c0;
      VY(1) =  c0;
      VZ(1) =  c1;

      VX(2) =  c0;
      VY(2) =  c2;
      VZ(2) =  c4;

      VX(3) =  c3;
      VY(3) =  c5;
      VZ(3) =  c4;

      VX(4) =  c6;
      VY(4) = -c7;
      VZ(4) =  c4;

      VX(5) = -c6;
      VY(5) = -c7;
      VZ(5) =  c4;

      VX(6) = -c3;
      VY(6) =  c5;
      VZ(6) =  c4;

      VX(7) = -c6;
      VY(7) =  c7;
      VZ(7) = -c4;

      VX(8) =  c6;
      VY(8) =  c7;
      VZ(8) = -c4;

      VX(9) =  c3;
      VY(9) = -c5;
      VZ(9) = -c4;

      VX(10) =  c0;
      VY(10) = -c2;
      VZ(10) = -c4;

      VX(11) = -c3;
      VY(11) = -c5;
      VZ(11) = -c4;

      VX(12) =  c0;
      VY(12) =  c0;
      VZ(12) = -c1;

%     for i=1:12
%      VX(i) = VX(1)/ru;
%      VY(i) = VX(1)/ru;
%      VZ(i) = VX(1)/ru;
%    end 

%------------------------
% define the vertex nodes
%------------------------

      x(1,1) = VX(1);  % first element
      y(1,1) = VY(1);
      z(1,1) = VZ(1);
      x(1,2) = VX(2);
      y(1,2) = VY(2);
      z(1,2) = VZ(2);
      x(1,3) = VX(3);
      y(1,3) = VY(3);
      z(1,3) = VZ(3);
%---
      x(2,1) = VX(1);
      y(2,1) = VY(1);
      z(2,1) = VZ(1);
      x(2,2) = VX(3);
      y(2,2) = VY(3);
      z(2,2) = VZ(3);
      x(2,3) = VX(4);
      y(2,3) = VY(4);
      z(2,3) = VZ(4);
%---
      x(3,1) = VX(1);
      y(3,1) = VY(1);
      z(3,1) = VZ(1);
      x(3,2) = VX(4);
      y(3,2) = VY(4);
      z(3,2) = VZ(4);
      x(3,3) = VX(5);
      y(3,3) = VY(5);
      z(3,3) = VZ(5);
%---
      x(4,1) = VX(1);
      y(4,1) = VY(1);
      z(4,1) = VZ(1);
      x(4,2) = VX(5);
      y(4,2) = VY(5);
      z(4,2) = VZ(5);
      x(4,3) = VX(6);
      y(4,3) = VY(6);
      z(4,3) = VZ(6);
%---
      x(5,1) = VX(1);
      y(5,1) = VY(1);
      z(5,1) = VZ(1);
      x(5,2) = VX(6);
      y(5,2) = VY(6);
      z(5,2) = VZ(6);
      x(5,3) = VX(2);
      y(5,3) = VY(2);
      z(5,3) = VZ(2);
%---
      x(6,1) = VX(2);
      y(6,1) = VY(2);
      z(6,1) = VZ(2);
      x(6,2) = VX(8);
      y(6,2) = VY(8);
      z(6,2) = VZ(8);
      x(6,3) = VX(3);
      y(6,3) = VY(3);
      z(6,3) = VZ(3);
%---
      x(7,1) = VX(3);
      y(7,1) = VY(3);
      z(7,1) = VZ(3);
      x(7,2) = VX(9);
      y(7,2) = VY(9);
      z(7,2) = VZ(9);
      x(7,3) = VX(4);
      y(7,3) = VY(4);
      z(7,3) = VZ(4);
%---
      x(8,1) = VX(4);
      y(8,1) = VY(4);
      z(8,1) = VZ(4);
      x(8,2) = VX(10);
      y(8,2) = VY(10);
      z(8,2) = VZ(10);
      x(8,3) = VX(5);
      y(8,3) = VY(5);
      z(8,3) = VZ(5);
%---
      x(9,1) = VX(5);
      y(9,1) = VY(5);
      z(9,1) = VZ(5);
      x(9,2) = VX(11);
      y(9,2) = VY(11);
      z(9,2) = VZ(11);
      x(9,3) = VX(6);
      y(9,3) = VY(6);
      z(9,3) = VZ(6);
%---
      x(10,1) = VX(6);
      y(10,1) = VY(6);
      z(10,1) = VZ(6);
      x(10,2) = VX(7);
      y(10,2) = VY(7);
      z(10,2) = VZ(7);
      x(10,3) = VX(2);
      y(10,3) = VY(2);
      z(10,3) = VZ(2);
%---
      x(11,1) = VX(2);
      y(11,1) = VY(2);
      z(11,1) = VZ(2);
      x(11,2) = VX(7);
      y(11,2) = VY(7);
      z(11,2) = VZ(7);
      x(11,3) = VX(8);
      y(11,3) = VY(8);
      z(11,3) = VZ(8);
%---
      x(12,1) = VX(3);
      y(12,1) = VY(3);
      z(12,1) = VZ(3);
      x(12,2) = VX(8);
      y(12,2) = VY(8);
      z(12,2) = VZ(8);
      x(12,3) = VX(9);
      y(12,3) = VY(9);
      z(12,3) = VZ(9);
%---
      x(13,1) = VX(4);
      y(13,1) = VY(4);
      z(13,1) = VZ(4);
      x(13,2) = VX(9);
      y(13,2) = VY(9);
      z(13,2) = VZ(9);
      x(13,3) = VX(10);
      y(13,3) = VY(10);
      z(13,3) = VZ(10);
%---
      x(14,1) = VX(5);
      y(14,1) = VY(5);
      z(14,1) = VZ(5);
      x(14,2) = VX(10);
      y(14,2) = VY(10);
      z(14,2) = VZ(10);
      x(14,3) = VX(11);
      y(14,3) = VY(11);
      z(14,3) = VZ(11);
%---
      x(15,1) = VX(6);
      y(15,1) = VY(6);
      z(15,1) = VZ(6);
      x(15,2) = VX(11);
      y(15,2) = VY(11);
      z(15,2) = VZ(11);
      x(15,3) = VX(7);
      y(15,3) = VY(7);
      z(15,3) = VZ(7);
%---
      x(16,1) = VX(7);
      y(16,1) = VY(7);
      z(16,1) = VZ(7);
      x(16,2) = VX(12);
      y(16,2) = VY(12);
      z(16,2) = VZ(12);
      x(16,3) = VX(8);
      y(16,3) = VY(8);
      z(16,3) = VZ(8);
%---
      x(17,1) = VX(8);
      y(17,1) = VY(8);
      z(17,1) = VZ(8);
      x(17,2) = VX(12);
      y(17,2) = VY(12);
      z(17,2) = VZ(12);
      x(17,3) = VX(9);
      y(17,3) = VY(9);
      z(17,3) = VZ(9);
%--- 
      x(18,1) = VX(9);
      y(18,1) = VY(9);
      z(18,1) = VZ(9);
      x(18,2) = VX(12);
      y(18,2) = VY(12);
      z(18,2) = VZ(12);
      x(18,3) = VX(10);
      y(18,3) = VY(10);
      z(18,3) = VZ(10);
%---
      x(19,1) = VX(10);
      y(19,1) = VY(10);
      z(19,1) = VZ(10);
      x(19,2) = VX(12);
      y(19,2) = VY(12);
      z(19,2) = VZ(12);
      x(19,3) = VX(11);
      y(19,3) = VY(11);
      z(19,3) = VZ(11);
%---
      x(20,1) = VX(11);
      y(20,1) = VY(11);
      z(20,1) = VZ(11);
      x(20,2) = VX(12);
      y(20,2) = VY(12);
      z(20,2) = VZ(12);
      x(20,3) = VX(7);
      y(20,3) = VY(7);
      z(20,3) = VZ(7);

%------------------------------------------
% compute the mid-points of the three sides
% of the 20 first-generation elements
%
% midpoints are numbered 4, 5, 6
%------------------------------------------

      Nelm = 20;

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
% compute the local element node coordinates
% for discretization levels 1 through Ndiv
%-------------------------------------------

      for i=1:Ndiv

       nm = 0;      % count the new elements arising by sub-division
                    % four element will be generated during each pass
       for j=1:Nelm  % over old elements

%---
% assign corner points to sub-elements
% these will become the "new" elements
%---

        nm = nm+1;

        xn(nm,1) = x(j,1);    %  first sub-element
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

        xn(nm,1) = x(j,4);      %  fourth sub-element
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

       end                      %  end of old-element loop

%--------------------------------------
% number of elements has increased
% by a factor of four
%--------------------------------------

       Nelm = 4*Nelm;

%---
% relabel the new points (xn -> x)
% and place them in the master list
%---

       for k=1:Nelm
        for l=1:6

         x(k,l) = xn(k,l);
         y(k,l) = yn(k,l);
         z(k,l) = zn(k,l);

%--- project onto the unit sphere

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

      end  %  of discretization-level loop

%-----------------------------------------
end
%-----------------------------------------


%-----------------------------------
% Generate a list of global nodes by looping 
% over all elementsand adding nodes not found
% in the current list.
%
% Fill in the connectivity table n(i,j) 
% containing node numbers of element points 1-6
%-----------------------------------

%------------------------------------
% three nodes of the first element are
% entered mannualy
%------------------------------------

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
% Iflag = 0 will signal a new global node
%---

      for i=2:Nelm        % loop over elements
       for j=1:3          % loop over element nodes

        Iflag = 0;

         for k=1:Npts
          if(abs(x(i,j)-p(k,1))<=eps) 
           if(abs(y(i,j)-p(k,2))<=eps)
            if(abs(z(i,j)-p(k,3))<=eps)

             Iflag  = 1;  % the node has been recorded previously
             n(i,j) = k;  % the jth local node of element i
                          % is the kth global node 
            end
           end
          end
         end
        
         if(Iflag==0) % record the node

          Npts = Npts+1;         % one more global node

          p(Npts,1) = x(i,j);
          p(Npts,2) = y(i,j);
          p(Npts,3) = z(i,j);

          n(i,j) = Npts;  % the jth local node of element i
                          % is the new global node 
         end 

       end 
      end          % of loop over elements


%----------------------------------
% Generate connectivity table: ne(i,j)
% for elements touching global node i
%
%  ne(i,j) ... ne(i,1) is the number of elements touching 
%                      the ith global node
%
%  ne(i,2:ne(i,1)+1) are the corresponding element labels
%----------------------------------

%---
% initialize
%---

      for i=1:Npts
       for j=1:4
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

         if(abs(p(i,1)-x(j,k)) <= eps)
          if(abs(p(i,2)-y(j,k)) <= eps)
           if(abs(p(i,3)-z(j,k)) <= eps)

            Icount = Icount+1;
            ne(i,1) = ne(i,1)+1;
            ne(i,Icount) = j;

           end
          end
         end

        end 
       end 
     
      end  %  of loop over global nodes


%------------------------------------------
% Generate connectivity table nbe(i,j)
%
%  nbe(i,j) .. label of element sharing 
%              the jth side of the ith element
%              for j = 1, 2, 3
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

      for i=1:Nelm       %  loop over elements
       jcount = 1;
       for j=4:6         %  loop over mid-points

        for k=1:Nelm          %  test element
         if(k~=i)             %  not a self-element
          for l=4:6             %  loop over mid-points

          if(abs(x(i,j)-x(k,l))<=eps)
           if(abs(y(i,j)-y(k,l))<=eps)
            if(abs(z(i,j)-z(k,l))<=eps)
             nbe(i,jcount) = k;
            end
           end
          end

          end
         end
        end %  of test element

        if(nbe(i,jcount)~=0)
         jcount = jcount+1;
        end

       end
      end        %  of loop over elements

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
% done
%-----

return

clear all

%====================================
%  Code edm_lob
%
%  Computation of the dimensionless
%  element difusion matrix for the
%  lobatto nodal base using
%  the lobatto quadrature
%==================================

%----------------------------
% lobatto interpolation nodes
%----------------------------

m = input('Enter the order of the polynomial expansion, m: ');

%---------------------------
% generate the spectral nodes
% and the integration weights
%---------------------------

xi(1) = -1.0;
w(1) = 2/(m*(m+1));
                                                                                
if(m >1)
  [tL, wL] =  lobatto(m-1);
  for i=2:m
   xi(i) = tL(i-1);
   w(i)  = wL(i-1);
  end
end

xi(m+1) = 1.0;

w(m+1)=w(1);

%----------------------------------------
% compute the node differentiation matrix
%----------------------------------------

for i=1:m+1
  for j=1:m+1

      %---
      if(i ~= j)
      %---
        ndm(i,j) = 1.0/(xi(j)-xi(i));
        for l=1:m+1
          if(l ~= j) ndm(i,j) = ndm(i,j)*(xi(j)-xi(l)); end
          if(l ~= i) ndm(i,j) = ndm(i,j)/(xi(i)-xi(l)); end
        end
      %---
      else
      %---
         ndm(i,i) = 0.0;
         for l=1:m+1
           if(i ~= l) ndm(i,i) = ndm(i,i) +1.0/(xi(i)-xi(l)); end
         end
      %---
      end
      %---
  end
end

%----------------------------
% compute the difusion matrix
% by the Lobatto quadrature
%----------------------------

for i=1:m+1    % loop over interpolation nodes
   for j=1:m+1        % loop over interpolation nodes

     dm(i,j) = 0.0;

     for l=1:m+1
       dm(i,j) = dm(i,j) + ndm(i,l)*ndm(j,l)*w(l);
     end

   end
end

%------
% print
%------

dm
det(dm)
inv(dm)

%-----
% done
%-----

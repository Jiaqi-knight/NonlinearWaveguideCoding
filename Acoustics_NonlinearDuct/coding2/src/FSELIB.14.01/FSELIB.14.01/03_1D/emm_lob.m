clear all

%===============================================
%  Code emm_lob
%
%  Computation of the dimensionless
%  element mass matrix
%  using the lobatto quadrature
%
%  m: order of polynomial expansion
%  k: order of Lobatto quadrature
%
%  setting m=k implements mass lumping
%
%================================================

%----------------------------
% lobatto interpolation nodes
%----------------------------

m = input('Enter the order of the polynomial expansion, m: ');

%---------------------------
% generate the spectral nodes
%---------------------------

xi(1) = -1.0;

if(m >1)
  [tL, wL] =  lobatto(m-1);
  for i=2:m
   xi(i) = tL(i-1);
  end
end

xi(m+1) = 1.0;

%----------------------------------------
% generate the lobatto integration points
% and corresponding weights
%----------------------------------------

disp ('Setting k=m implements mass lumping')
k = input('Enter the order of the lobatto quadrature, k: ');

z(1) = -1.0; z(k+1)=1.0;
w(1) = 2/(k*(k+1)); w(k+1)=w(1);

if(k>1)
  [zL, wL] =  lobatto(k-1);
  for i=2:k
   z(i) = zL(i-1); w(i) = wL(i-1);
  end
end

%--------------------------
% compute the mass matrix
% by the Lobatto quadrature
%--------------------------

for ind=1:m+1    % loop over interpolation nodes
   for jnd=1:m+1        % loop over interpolation nodes

  % Generate the integrand of the mass matrix

%--------
      for j=1:k+1
%--------

       f(j) = 1.0;

       for i=1:m+1
        if(i ~= ind)
         f(j) = f(j)*(z(j)-xi(i))/(xi(ind)-xi(i));
        end
        if(i ~= jnd)
         f(j) = f(j)*(z(j)-xi(i))/(xi(jnd)-xi(i));
        end
       end

%--------
      end  % of loop over j
%--------

  % perform the quadrature

     msmat = 0.0;

     for i=1:k+1
       msmat = msmat + f(i)*w(i);
     end

    mm(ind,jnd) = msmat;

   end
end

%------
% print
%------

mm

%-----
% done
%-----

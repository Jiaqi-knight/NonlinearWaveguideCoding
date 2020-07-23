function [x, ...
          l,u,det, ...
          Istop ] = gel (n,a,rhs,Iwlpvt,Isym)

%============================================
% Solution of the nxn linear system:
%
%    a x = rhs
%
% by Gauss elimination with row pivoting
%
% If Iwlpvt = 1, row pivoting is enabled
% If Isym   = 1, matrix a is symmetric
%
% c: extended coefficient matrix
% l: lower triangular matrix
% u: upper triangular matrix
%
% det: determinant of a
%============================================

%-----------
% initialize
%-----------

Istop  = 0;          % error flag
Icount = 0;          % counts row interchanges in pivoting
eps = 0.00000000001; % tolerance

%-----------
% pivoting is not done for a symmetric system
%-----------

if(Isym==1)
 disp(' gel: system is symmetric; pivoting is disabled')
 Iwlpvt = 0;
end

%--------
% prepare
%--------
                                                                                
na = n-1; n1 = n+1;

%-------------------
% Initialize l and c
%-------------------

for i=1:n
  for j=1:n
    l(i,j) = 0.0; c(i,j) = a(i,j);
  end
  c(i,n1) = rhs(i);
end
                                                                                
%---------------------
% Begin row reductions
%---------------------

for m=1:na   % outer loop for working row

   ma = m-1; m1 = m+1;

%-------------------------
% Pivoting module:
%
% search the ith column
% for the largest element
%-------------------------

if(Iwlpvt==1)   % pivoting will be done if Iwlpvt=1
                                                                                
  Ipv = m; pivot = abs(c(m,m));

    for j=m1:n
      if(abs(c(j,m)) > pivot) 
       Ipv = j; pivot = abs(c(j,m));
      end
    end

    if(pivot < eps) 
      disp ('gel: trouble in station 1')
      Istop = 1;
      return
    end

%--------------------------------------
% Switch the working row with
% the row containing the pivot element
% Also switch rows in l
%--------------------------------------

    if(Ipv ~= m) 
                                                                                
     for j=m:n1
       save = c(m,j); c(m,j) = c(Ipv,j); c(Ipv,j) = save;
     end
     for j=1:ma
       save = l(m,j); l(m,j) = l(Ipv,j); l(Ipv,j) = save;
     end

     Icount = Icount+1;    % increase the pivoting counter
                                                                                
   end

end   % end of pivoting module for Iwlpvt
                                                                                
%------------------------------------------
% reduce column i underneath element c(m,m)
%------------------------------------------
                                                                                
 for i=m1:n
                                                                                
  l(i,m) = c(i,m)/c(m,m);

  if(Isym==1) 
   l(i,m) = c(m,i)/c(m,m); ilow = i; 
  else
   l(i,m) = c(i,m)/c(m,m); ilow = m1; 
  end

  c(i,m) = 0.0;

  for j=ilow:n1
    c(i,j) = c(i,j)-l(i,m)*c(m,j);
  end

 end

%---
% end of outer loop for working row:
%---

end
                                                                                
%---------------------------------
% check the last diagonal element
% for singularity
%--------------------------------
                                                                                
if(abs(c(n,n)) < eps)

   disp('gel: trouble in station 2')
   Istop = 1;
   return;

end
                                                                                
%----------------------
% complete the matrix l
%----------------------

for i=1:n
  l(i,i)=1.0;
end

%--------------------
% define the matrix u
%--------------------
                                                                                
for i=1:n
  for j=1:n
    u(i,j) = c(i,j);
  end
end

%-----------------------------------
% perform back-substitution to solve
% the reduced system
% using the upper triangular matrix c
%------------------------------------
                                                                                
x(n) = c(n,n1)/c(n,n);
                                                                                
for i=na:-1:1
  sum=c(i,n1);
  for j=i+1:n
    sum=sum-c(i,j)*x(j);
  end
  x(i)=sum/c(i,i);
end
                                                                                
%-----------------------
% compute the determinant as:
%
% det(a) = (+-) det(l)*det(u)
%-----------------------
                                                                                
det = 1.0;
                                                                                
for i=1:n
   det=det*c(i,i);
end
                                                                                
if(Iwlpvt == 1)
  for i=1:Icount
    det = -det;
  end
end

%-----
% done
%-----

return

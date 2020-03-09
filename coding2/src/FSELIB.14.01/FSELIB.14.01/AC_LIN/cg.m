function [sln] = cg (n, a, rhs)

%-------------------------------------------------------
%  Solution of a linear symmetric positive-definite
%  system by the method of conjugate gradients
%  The search vectors are chosen by the method of
%  Hestenes and Steifel (1952)
%
%  SYMBOLS:
%
%  a .... symmetric positive definite matrix
%  n .... size (rows/columns) of matrix a
%  rhs .. right hand side vector (e.g. b, as in Ax=b)
%  x .... evolving solution vector
%  sln .. final solution vector
%  p .... search directions
%  r .... residual vectors
%  alpha. scale parameter for solution update
%  beta.. scale parameter for search direction
%
%-------------------------------------------------------
  
%---------------------------------------------------
% set the initial values of the vectors x and r (step 0)
% set the first-step value of p
%---------------------------------------------------
                                                     
for i=1:n
  x0(i) = 0.0; r0(i) = rhs(i); p(1,i) = rhs(i);
end

%----------------------------
% form the sums used in alpha
%----------------------------
                              
alpha_num = 0.0; alpha_den = 0.0;

for i=1:n
  alpha_num = alpha_num + r0(i)*r0(i);
  for j=1:n
    alpha_den = alpha_den + p(1,i)*a(i,j)*p(1,j);
  end
end

alpha = alpha_num/alpha_den;

%---------------------------------------------------
% set first step values of alpha and vectors x and r
%---------------------------------------------------
                                                     
for i=1:n
  x(1,i) = x0(i)+alpha*p(1,i);
  r(1,i) = r0(i);
  for j=1:n
    r(1,i) = r(1,i)-alpha*a(i,j)*p(1,j);
  end
end

%-------------------------------------------
% loop through the remaining search vectors
% 2 to n, and compute
% alpha, beta, and vectors p, x, and r
%-------------------------------------------
                                             
for k=2:n       %  outer loop over search directions
                                                             
%---
% sums used in beta
%---
                    
beta_num = 0.0; beta_den = 0.0D0;
                       
for i=1:n
  beta_num = beta_num + r(k-1,i)^2;
  if(k==2)
    beta_den = beta_den + r0(i)^2;
  else
    beta_den = beta_den + r(k-2,i)^2;
  end
end

beta = beta_num/beta_den;
                               
for i=1:n
  p(k,i) = r(k-1,i)+beta*p(k-1,i);
end

%---
% compute the sums used in alpha
%---
                                 
alpha_num = beta_num; alpha_den = 0.0;
                           
for i=1:n
  for j=1:n
   alpha_den = alpha_den + p(k,i)*a(i,j)*p(k,j);
 end
end

alpha = alpha_num/alpha_den;
                                  
%---
% compute the k'th iterations of vectors x and r
%---
                                                 
for i=1:n
   x(k,i) = x(k-1,i)+alpha*p(k,i);
   r(k,i) = r(k-1,i);
   for j=1:n
      r(k,i)=r(k,i)-alpha*a(i,j)*p(k,j);
   end
end
               
end              % end of outer loop

%---------------------
% Extract the solution
%---------------------
                       
for i=1:n
  sln(i) = x(n,i);
end
                        
return

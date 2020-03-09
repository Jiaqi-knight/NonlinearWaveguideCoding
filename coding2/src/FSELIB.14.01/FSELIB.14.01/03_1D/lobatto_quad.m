function integral = int_lob(k)

%====================================
% int_lob
%
% Numerical integration of a function
% in the interval [-1, 1]
% using the Lobatto quadrature
%====================================

% k = input('Please enter the order of the quadrature, k: ');

%------------------------------------------------
% make up the base points and integration weights
%------------------------------------------------

z(1) = -1.0;
w(1) = 2.0/(k*(k+1));

if(k>1) 

  [zL, wL]= lobatto(k-1);

  for i=2:k
   z(i) = zL(i-1);
   w(i) = wL(i-1);
  end

end

z(k+1)=1.0;
w(k+1)=w(1);

%-----------------------------
% generate the function values
%-----------------------------

for i=1:k+1
  xi = z(i);
  f(i) = xi^12;
  f(i) = xi^11;  % examples
  f(i) = xi^10;
end

%-----------------------
% perform the quadrature
%-----------------------

integral = 0.0;

for i=1:k+1
  integral = integral + f(i)*w(i);
end

%-----
% done
%-----

return

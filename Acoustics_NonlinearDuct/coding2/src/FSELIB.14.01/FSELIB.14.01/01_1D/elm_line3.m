function xe = elm_line3 (x1,x2,n,ratio);

%===========================================
% FSELIB
%
% Symmetric discretization of a line segment
% into a graded mesh of n elements
% subtended between the left end-point x1
% and the right end-point x2
%
% n is odd: 1,3,5,...
%
% The variable "ratio" is the ratio
% of the size of the mid-elements
% to the size of the first element
%============================================

%----------
% 1 element
%----------

if(n==1)
  xe(1) = x1;
  xe(2) = x2;
  return;
end

%--------------
% many elements
%--------------

if(ratio==1)
  alpha  = 1.0;
  factor = 1.0/n;
else
  texp   = 2.0/(n-1.0);
  alpha  = ratio^texp;
  tmp1   = (n+1.0)/2.0;
  tmp2   = (n-1.0)/2.0;
  factor = (1.0-alpha)/(2.0-alpha^tmp1-alpha^tmp2);
end

deltax = (x2-x1) * factor;   % length of first element

%-----------
% discretize
%-----------

xe(1) = x1;    % first point

for i=2:(n+3)/2
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax*alpha;
end

deltax = deltax/alpha^2;

for i=(n+5)/2:n+1
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax/alpha;
end

%-----
% done
%-----

return;

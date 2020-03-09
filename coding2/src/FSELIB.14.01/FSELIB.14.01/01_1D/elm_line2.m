function xe = elm_line2(x1,x2,n,ratio);

%===========================================
% FSELIB
%
% Symmetric discretization of a line segment
% into a graded mesh of n elements
% subtended between the left end-point x1
% and the right end-point x2
%
% n is 1 or even: 2,4,6,...
%
% The variable "ratio" is the ratio
% of the size of the mid
% elements to the size of the first element
%===========================================

%------------
% one element
%------------

if(n==1)
  xe(1) = x1;
  xe(2) = x2;
  return;
end

%-------------
% two elements
%-------------

if(n==2)
  xe(1) = x1;
  xe(2) = 0.5*(x1+x2);
  xe(3) = x2;
  return;
end

%--------------
% many elements
%--------------

xe = zeros(n+1,1);
nh = n/2;
xh = 0.5*(x1+x2);   % mid-point

if(ratio==1)
  alpha  = 1.0;
  factor = 1.0/nh;
else
  texp   = 1.0/(nh-1.0);
  alpha  = ratio^texp;
  factor = (1.0-alpha)/(1.0-alpha^nh);
end

deltax = (xh-x1) * factor;   % length of first element

%-----------
% discretize
%-----------

xe(1) = x1;    % first point

for i=2:nh+1
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax*alpha;
end

deltax = deltax/alpha;

for i=nh+2:n+1;
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax/alpha;
end

%-----
% done
%-----

return;

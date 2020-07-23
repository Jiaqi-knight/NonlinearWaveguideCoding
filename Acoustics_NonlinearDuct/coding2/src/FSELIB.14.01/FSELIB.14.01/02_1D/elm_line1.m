function xe = elm_line1 (x1,x2,n,ratio);

%===============================================
% FSELIB
%
% Unsymmetrical discretization of a line segment
% into a graded mesh of n elements
% subtended between the left end-point x1
% and the right end-point x2
%
% The variable "ratio" is the ratio
% of the length of the last
% element to the length of the first element.
%
% xe: element end-nodes
%===============================================

%------------
% one element
%------------

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
  texp = 1/(n-1);
  alpha = ratio^texp;
  factor = (1.0-alpha)/(1.0-alpha^n);
end

deltax=(x2-x1)*factor;   % length of the first element

xe(1)=x1;    % first point

for i=2:n+1
  xe(i)  = xe(i-1)+deltax;
  deltax = deltax*alpha;
end

%-----
% done
%-----

return;

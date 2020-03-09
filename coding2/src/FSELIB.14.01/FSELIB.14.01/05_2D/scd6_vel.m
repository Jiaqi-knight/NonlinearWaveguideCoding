function [u, v] = scd6_vel(U,a,x,y)

%===================
% evaluation of nodal velocity
%===================

  as = a^2;
  rs = x^2+y^2;
  rq = rs^2;

  u = U*(1.0+as/rs - 2.0*as*x^2/rq);
  v = U*( -2.0*as*x*y/rq);

%---
% done
%---

return;

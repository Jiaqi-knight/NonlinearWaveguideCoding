function x = thomas (n,a,b,c,s)

%==================================================
% FSELIB
%
% Thomas algorithm for a tridiagonal system
%
% n:     system size
% a,b,c: diagonal, superdiagonal,
%        and subdiagonal elements
% s:     right-hand side
%==================================================

%--------
% prepare
%--------

na = n-1;

%------------------------------
% reduction to upper bidiagonal
%------------------------------

d(1) = b(1)/a(1);
y(1) = s(1)/a(1);

for i=1:na
  i1 = i+1;
   den   = a(i1)-c(i1)*d(i);
   d(i1) = b(i1)/den;
   y(i1) = (s(i1)-c(i1)*y(i))/den;
end

%------------------
% back substitution
%------------------

x(n) = y(n);

for i=na:-1:1
  x(i) = y(i)-d(i)*x(i+1);
end

%-----
% done
%-----

return;

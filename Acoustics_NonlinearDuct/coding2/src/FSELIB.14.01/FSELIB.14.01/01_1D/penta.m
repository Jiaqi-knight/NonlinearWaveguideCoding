function x = penta(n,a,b,c,d,e,s)

%==================================================
% FSELIB
%
% modified Thomas algorithm for
% a pentadiagonal linear system
%
% n:       system size
% a,b,c,d,e: diagonal, superdiagonal, super-super,
%        subdiagonal, and sub-sub  elements
% s:     right-hand side
%==================================================

na = n-1;
nb = n-2;

%-----------------------------
% reduction to quadra-diagonal
%-----------------------------

for i=1:nb
  c1(i) = c(i);
end

c1(na)=0.0; c1(n)=0.0;

a1(1) = a(1); b1(1) = b(1); s1(1) = s(1);

d1(2) = d(2); a1(2) = a(2); b1(2) = b(2); s1(2) = s(2);

for i=2:na
  i1 = i+1;
  w = e(i1)/d1(i);
  a1(i1) = a(i1) - w*b1(i);
  b1(i1) = b(i1) - w*c1(i);
  d1(i1) = d(i1) - w*a1(i);
  s1(i1) = s(i1) - w*s1(i);
end

%-----------------------------
% Generate the system T x= s2
%-----------------------------
                                                                                
b2(1) = b1(1)/a1(1);
c2(1) = c1(1)/a1(1);
s2(1) = s1(1)/a1(1);
                                                                                
for i=1:na
  i1 = i+1;
  den    =  a1(i1)-d1(i1)*b2(i);
  b2(i1) = (b1(i1)-d1(i1)*c2(i))/den;
  c2(i1) =  c1(i1)/den;
  s2(i1) = (s1(i1)-d1(i1)*s2(i))/den;
end

%------------------
% back substitution
%------------------
                                                                                
x(n) = s2(n);

x(n-1)= s2(n-1)-b2(n-1)*x(n);
                                                                                
for i=nb:-1:1
   x(i) = s2(i)-b2(i)*x(i+1)-c2(i)*x(i+2);
end

%-----
% done
%-----

return;

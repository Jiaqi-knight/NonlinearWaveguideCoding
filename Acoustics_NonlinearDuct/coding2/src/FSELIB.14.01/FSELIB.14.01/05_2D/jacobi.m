function jac = jacobi(a,b,n,t)

%=========================================
% evaluates by recursion Jacobi polynomial
%
% J^(a,b)_n
%=========================================

jac=1.0;

if(n>0)

al = 0.5*(a+b+2); be= (b-a)/(a+b+2);
p(1)= al*(t-be);

if(n>1)

  al = (2+a+b+1)*(2+a+b+2)/(4*(a+b+2) );
  be = (b^2-a^2)/((2+a+b)*(2+a+b+2));
  ga = (1+a)*(1+b)*(2+a+b+2)/(2*(a+b+2)*(2+a+b));

  p(2)= al*(t-be)*p(1) - ga;

if(n>2)

  for i=2:n-1
    al = (2*i+a+b+1)*(2*i+a+b+2)/(2*(i+1)*(i+a+b+1) );
    be = (b^2-a^2)/((2*i+a+b)*(2*i+a+b+2));
    ga = (i+a)*(i+b)*(2*i+a+b+2)/((i+1)*(i+a+b+1)*(2*i+a+b));
    p(i+1)= al*(t-be)*p(i) - ga*p(i-1);
  end

end
end
  jac=p(n);
end

%=== Lobatto
% Lo = 0.5*(n+2)*jac
%===
% Lo1 = 3*t
% Lo5 = 1/8  * (693*t^4-630*t^2+105)*t
% Lo6 = 1/16 * (3003*t^6-3465*t^4+945*t^2-35)
%
%=== Legendre
% L = jac
% L6 =  1/16 * (231*t^6 - 315*t^4 + 105*t^2 - 5)
%=====

%-----
% done
%-----

return

%laguerre
% Computes the associated laguerre polynomials y(x) = laguerre(x,n,L).
% n is an integer 1,2,3,4... and L is 0,1,2....(n-1).
%
% LAST UPDATED by Andy Frennch Dec 2011

function y = laguerre(x,n,L)

y = zeros(size(x));
for k=0:n-L-1
    a1 = factorial(L+n);
    a2 = factorial(2*L+1+k);
    a3 = factorial(n-L-1-k);
    a4 = factorial(k);
    y = y + a1*((-x).^k )/( a2*a3*a4 );
end

%

function test_laguerre

x = linspace( -20*pi,20*pi, 100 );
n = 7;
L = 1;
y = laguerre(x,n,L);
plot(x,y);
xlabel('x')
ylabel('y')
title(['laguerre(x,',num2str(n),',',num2str(L),')']);

%End of code
clc
clear

A = magic(3)

syms x

f(x) = x;

B = funm(A,f)

double(B)
% inline(B)
% 
syms x;y=besselj(1,x);y1=inline(y,x);
Y = 

Y(1,1)
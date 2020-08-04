
%% written by saeid zeynali- saedznl1374@gmail.com - msc student of power systems engineering at university of Tabriz - Iran

% this is a matlab function that calculates the lagrange multipliers of
% the maximum ntropy problem (example provided)
% the enropy function is considered to be in form of p(x)=(lambda0+lambda1*x+lambda2*x.^2+lamda3*x^.3+lamda4*x.^4+lamda5*x^5+lamda6*x^6);
% which is ploynominal so the probibility density function would be in form of phi(x)=exp(p(x))
% with known lamda values obtained from this function

%% using instructions 
%xmin is minimum value of your probability density ,and xmax is the maximum
%of it. varargin (variable input argument) is the moments of your variables
%wich you want to find the problity density of it ,varargin must be at least
%two and maximum six. i'd like to give an exaxmple;

%% example: we have the first four non-central moments of a density function as m1,m2,m3,m4
%% and we want remake the function by maximum antripy method.
clc
clear
%fist write a function handle like --> 
xmin=0;xmax=1;
[F,flag,phi]=Maximum_Entropy(xmin,xmax,6,4,2,1) %varargin are the moments
f=@(x) phi(x,F);
% f is the probibility density function ;

%% plot
x=xmin:0.001:xmax;
y=f(x);
plot(x,y);
%% cautions 
% if you have 6 moments you should write phi=@(x,la) exp(-la(1)-la(2)*x-la(3)*x.^2-la(4)*x.^3-la(5)*x.^4-la(6)*x.^5-la(7)*x.^6);
% if problem is not converging increase nfe unde each case inside the
% function
% if its too slow decrese nfe under each case

%ME2
% This script shows how to use the function ME_DENS2
% in the case of the quartic distribution. (see Example 2.)
clc
clear
xmin=-1; xmax=1; dx=0.01; % define the x axis
x=[xmin:dx:xmax]';
mu=[0.1,.3,0.1,.15]'; % define the mu values
[lambda,p,entr]=reni2(mu,x,1);
disp([mu;lambda;entr(length(entr))]')
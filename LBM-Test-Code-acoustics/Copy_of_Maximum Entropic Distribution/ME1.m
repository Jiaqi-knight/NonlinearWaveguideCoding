%
% Permission was provided by the author, Ali Mahammad-Djafari, 
% to modify and distribute this code with the informME 
% package. Contact the original author directly for use outside 
% this package.
%
% Author's website:
% http://djafari.free.fr/index.htm
% 
% Author's paper with source code:
% Mohammad-Djafari A. (1992) A Matlab Program to Calculate the Maximum Entropy Distributions. In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds) Maximum Entropy and Bayesian Methods. Fundamental Theories of Physics (An International Book Series on The Fundamental Theories of Physics: Their Clarification, Development and Application), vol 50. Springer, Dordrecht
%

%WILL BE MODIFIED BY JIAQI WANG, FROM SJTU
%----------------------------------
%ME1
% This script shows how to use the function ME_DENS1
% in the case of the Gamma distribution. (see Example 1.)
clc;clear;close all
xmin=0.0001; xmax=1; dx=0.01; % define the x axis
x=[xmin:dx:xmax]';
mu=[0.3,-1.5]'; % define the mu values
[lambda,p,entr]=me_dens1(mu,x);
alpha=-lambda(3); beta=lambda(2);
m=(1+alpha)/beta; sigma=m/beta;
disp([mu' alpha beta m sigma entr(length(entr))])
%----------------------------------

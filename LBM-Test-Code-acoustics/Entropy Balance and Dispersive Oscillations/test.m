clc
clear
close all
%LBM Lattice Boltzmann Method
% LBM(MC,V,NX,TM,MOV,ESPIL,LIMIT,NORM) solves the shock tube problem
% with number of lattice
% sites NX with separation of one and viscosity V upto TM time steps of
% length one using method
% MC:
% 0 − LBGK polynomial equilibria,
% 1 − LBGK entropic equilibria,
% 2 − ELBM Newton iterations,
% 3 − ELBM Bisection method
% and entropic limiter
% LIMIT:
% 0 − No limiter
% 1 − Median filtering
% Output of a movie of the simulation can be controlled with
% MOV:
% 0 − No movie,
% 1 − With movie
% The accuracy of the root finding for ELBM can be controlled with
% EPSIL in all cases the root found will result in entropy production.
% Convergence is partly based on | | feq − f|| , the choice of the norm is
% controlled with
% NORM:
% 0 − L1 norm
% 1 − Entropic Norm;
V=0.5;%beta=1 %粘性
LBM(1,V,1000,20000,1,0.0000001,0,1)


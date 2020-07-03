%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 26th, 2014
%    Last update: May 26th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

clear all
close all
clc

Q = 5;

allvelocities = [0;1;-1;3;-3];
velocities = [0;1;3];

%allvelocities = [0;1;-1;2;-2;3;-3;5;-5];
%velocities = [0;1;2;3;5];

T0 = T0closurerelation1D(Q,velocities);

W = computeWi1D(Q,allvelocities,T0(2));

[scheme1D,cs,cssq,T0]=initializeELBMscheme(1,5);

D = 3;

[schemeQ125,QD125] = constructentropicDnscheme(D,scheme1D);

pruneflag = 1;

prune = [2;9;10;11;18;19;27];
[schemeQ15,QD15] = constructentropicDnscheme(D,scheme1D,pruneflag,prune);

prune = [3;9;10;11;18;19;27];
[schemeQ19,QD19] = constructentropicDnscheme(D,scheme1D,pruneflag,prune);

prune = [9;10;11;18;19;27];
[schemeQ27,QD27] = constructentropicDnscheme(D,scheme1D,pruneflag,prune);

prune = [10;11;18;19];
[schemeQ41,QD41] = constructentropicDnscheme(D,scheme1D,pruneflag,prune);

D = 2;

[schemeQ25,QD25] = constructentropicDnscheme(D,scheme1D);

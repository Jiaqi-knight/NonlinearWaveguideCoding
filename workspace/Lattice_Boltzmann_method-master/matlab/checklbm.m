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
[schemeQ41,cs,cssq,T0]=initializeELBMscheme(3,41);

figure;
for k=1:41
quiver3(0,0,0,schemeQ41(k,1),schemeQ41(k,1),schemeQ41(k,1),'.'); hold on;
% quiver3(...
% 	downsample(xyz(:,:,k,7),mxds),...
% 	downsample(xyz(:,:,k,8),mxds),...
% 	downsample(xyz(:,:,k,9),mxds),...
% 	downsample(cobase(:,:,k,7), mxds),...
% 	downsample(cobase(:,:,k,8), mxds),...
% 	downsample(cobase(:,:,k,9), mxds));
end
title 'Nhat';
% Create Nhat.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 3.0 3.0])
%print -depsc figures/Nhat.eps



D = 2;

[schemeQ25,QD25] = constructentropicDnscheme(D,scheme1D);

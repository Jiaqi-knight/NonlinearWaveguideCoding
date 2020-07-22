%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme
%D1Q7
%Langrange multipliers
%ELBM
%% Initialisation
clc
clear

% % Assumes D1Q7 and 7th trajectory is rest particle
% clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
c_scale = [-0.5 1.5];        % Velocity scale for plotting

% Time
T = 200;        % Number of time steps

% Lattice dimensions
NXMAX = 1000;  % Number of x lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:NXMAX;

% Discrete velocities for D1Q5
D=1;Q=7
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
% Initialise matrices (N x N (x 9))
u = zeros(NXMAX,1);       % Velocity field
rho = ones(NXMAX,1);       % Density field
k=2*pi/50;x=1:50;
% rho(x)=1+0.01*(1-cos(k*x));
% u(x)=-0.01/sqrt(3)*(1-cos(k*x));
% rho(1:NXMAX/2)=3;rho(NXMAX/2+1:end)=1;

[feq] = entropyEquilibrium(NXMAX,D,Q,T0,rho,scheme(:,2).',scheme(:,1).',u);
f=feq;%接下来是为了迭代初值

omega=2.1;

x=1:50;
rho(1:50)=1+0.01*(1-cos(k*x));
                          

figure
%% LBM loop
lambda=100;
W=2*pi/lambda;
omega=1;
tic
for t = 1 : 20
    plot(rho);     
    pause(0.01)
    % Compute collision
    [f] = Collide(NXMAX,D,Q,T0,rho,scheme(:,2),scheme(:,1),u,cs,feq,NXMAX,f,omega);
    % Stream populations
    f = Stream(Q,f,NXMAX,scheme(:,1));
    % Find macroscopic quantities
    [rho,u] = Macroscopic(f,scheme(:,1));
    % Call visualisation
    u(:)=0.5; %模拟均匀流条件；%需要小于0.8

end
toc


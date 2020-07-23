%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme
%D1Q7
%Langrange multipliers
%ELBM

%Add:
%changable local temperture

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
Time = 2000;        % Number of time steps

% Lattice dimensions
NXMAX = 300;  % Number of x lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:NXMAX;

% Discrete velocities for D1Q5
D=1;Q=7;
% Initialise matrices (N x N (x 9))
u = zeros(NXMAX,1);       % Velocity field
rho = ones(NXMAX,1);       % Density field
k=2*pi/50;x=1:50;
% rho(x)=1+0.01*(1-cos(k*x));
% u(x)=-0.01/sqrt(3)*(1-cos(k*x));
% rho(1:NXMAX/2)=3;rho(NXMAX/2+1:end)=1;
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);


[feq] = entropyEquilibrium(NXMAX,D,Q,T0,rho,scheme(:,2).',scheme(:,1).',u); %
f=feq;%接下来是为了迭代初值

omega=1.9;

% x=1:50;
% rho(1:50)=1+0.01*(1-cos(k*x));
                          

figure
%% LBM loop
% lambda=100;
% W=2*pi/lambda;
% omega=1;
tic
for t = 1 : 100
    plot(rho);     
    pause(0.01)
    % Compute collision
    [f] = Collide(NXMAX,D,Q,T0,rho,scheme(:,2),scheme(:,1),u,cs,feq,NXMAX,f,omega);
    % Stream populations
    f = Stream(Q,f,NXMAX,scheme(:,1));
    % Find macroscopic quantities
    [rho,u] = Macroscopic(f,scheme(:,1));
    % Call visualisation
     rho(NXMAX/2)=1+0.01*(cos(k*t));
%     u(:)=0; %模拟均匀流条件；%需要小于0.8
end
toc


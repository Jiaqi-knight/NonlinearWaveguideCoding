%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme

%% Initialisation
% clc
% % Assumes D1Q5 and 5th trajectory is rest particle
% clear
set(0,'DefaultFigureWindowStyle','docked')
close all


N = 300;  % Number of x lattice sites
omega = 1.9;    % Relaxation frequency
   
% Discrete velocities for D1Q3
D=1;Q=3
[scheme,cs]=initializeELBMscheme(D,Q);

f=zeros(N,Q);                                 
feq=zeros(N,Q);
f(:,:)=1/3;   
u = zeros(N, 1);

% rho_l = 0.01;   % initial disturbance
% f(N/2,3) = rho_l;

% rho(N/2)=0.01
  k=2*pi/50;
%x=1:50;
% rho(x)=1+0.01*(1-cos(k*x));
% u(x)=-0.01/sqrt(3)*(1-cos(k*x));
% rho(1:NXMAX/2)=3;rho(NXMAX/2+1:end)=1;

% [feq] = Equilibrium(rho,scheme(:,2).',scheme(:,1).',u);
% f=feq;%接下来是为了迭代初值



                          

figure
%% LBM loop
% lambda=10;
% W=2*pi/lambda;
tic
for t = 1 : 100
  
    
    f = Stream(f,N,scheme(:,1));
    [rho,u] = Macroscopic(f,scheme(:,1));
    
    rho(N/2)=1+0.01*(cos(k*t));
    % Compute collision
    [f] = Collide(rho,scheme(:,2),scheme(:,1),u,f,omega);
    % Stream populations
    % Find macroscopic quantities
    % Call visualisation
    %
%     u(:)=0; %模拟均匀流条件；%需要小于0.8
    plot(rho);     
    pause(0.01)  
end
toc


%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme
%D1Q7
%Langrange multipliers
%ELBM

%Add:
%changable local temperture

%
clc
clear

% % Assumes D1Q7 and 7th trajectory is rest particle
% clear
set(0,'DefaultFigureWindowStyle','docked')
% close all


% Lattice dimensions
NX = 300;  % Number of x lattice sites
omega = 1;    % Relaxation frequency


% Discrete velocities for D1Q7
D=1;Q=7;
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
% Initialisation
f=zeros(NX,Q);                                 
feq=zeros(NX,Q);
f(:,:)=1/Q;   
u = zeros(NX, 1);


k=2*pi/50;
% rho(x)=1+0.01*(1-cos(k*x));
% u(x)=-0.01/sqrt(3)*(1-cos(k*x));
% rho(1:NXMAX/2)=3;rho(NXMAX/2+1:end)=1;


% [feq] = entropyEquilibrium(NX,D,Q,T0,rho,scheme(:,2).',scheme(:,1).',u); %
% f=feq;%接下来是为了迭代初值


                          

figure
%% LBM loop

tic
for t = 1 : 100
    f = Stream(Q,f,NX,scheme(:,1));
    [rho,u] = Macroscopic(f,scheme(:,1));
%     u=u+0.011; %模拟均匀流条件；%需要小于0.8

%     rho(NX/2)=1+0.01*(cos(k*t));
    [f] = Collide(NX,D,Q,T0,rho,scheme(:,2),scheme(:,1),u,cs,feq,NX,f,omega);
   
    plot(rho);     
    pause(0.01)
end
toc


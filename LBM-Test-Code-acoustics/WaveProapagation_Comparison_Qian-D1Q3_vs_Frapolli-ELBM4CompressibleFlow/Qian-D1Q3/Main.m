%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme

%% Initialisation
clc
clear
% close all


NX = 1000;  % Number of x lattice sites
omega = 2;    % Relaxation frequency
   
% Discrete velocities for D1Q3
D=1;Q=3;
[scheme,cs]=initializeELBMscheme(D,Q);
% Initialisation
f=zeros(NX,Q);                                 
feq=zeros(NX,Q);
f(:,:)=1/Q;   
u = zeros(NX, 1);

% rho_l = 0.01;   % initial disturbance
% f(N/2,3) = rho_l;
f(NX/2,1)=0.01

k=2*pi/50;


figure
%% LBM loop
tic
for t = 1 : 1000
    f = Stream(f,NX,scheme(:,1));   
%     k=2*pi/50;
%     f(NX/2,1)=f(NX/2,1)+0.1*(cos(k*t));
    [rho,u] = Macroscopic(f,scheme(:,1));
    %u=u+0.001; %模拟均匀流条件；%需要小于0.8

%     rho(NX/2)=1+0.01*(cos(k*t));
    % Compute collision
    [f] = Collide(rho,scheme(:,2),scheme(:,1),u,f,omega);
    plot(rho);     
    pause(0.01)  
end
toc


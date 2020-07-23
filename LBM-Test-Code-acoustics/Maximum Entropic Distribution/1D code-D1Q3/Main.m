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

% Visuals
c_scale = [-0.5 1.5];        % Velocity scale for plotting

% Time
T = 10000;        % Number of time steps

% Lattice dimensions
NXMAX = 1000;  % Number of x lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:NXMAX;

% Discrete velocities for D1Q5
D=1;Q=3

%% 2 - Set physical parameters (macro)
physical_sound_velocity = 340; % [m/s]
physical_density = 1.2; % [kg/m^3]
physical_dimension_max_x = 0.3; % [m]
physical_dimension_max_y = 0.072; % [m]
% voxel is a term to express a volume decribed in a pixel: volume + pixel = voxel
dimension_x_voxel = physical_dimension_max_x/number_columns_lattice; % defining dimension x in voxel
lattice_time_step = (1/sqrt(3))*dimension_x_voxel/physical_sound_velocity;

%% 3 - Set lattice parameters (meso - lattice unities)
frequency_relaxation = 1.9; % to 1.5e-5 physcosity 1.9998; 860e-5 = 1.9
time_relaxation = 1/frequency_relaxation;
lattice_average_density = 1;
lattice_sound_speed = 1/sqrt(3);
lattice_sound_speed_pow_2 = lattice_sound_speed^2;
lattice_viscosity = lattice_sound_speed_pow_2*(1/frequency_relaxation-0.5);
physical_viscosity = lattice_viscosity*(dimension_x_voxel^2)/lattice_time_step; % [m^2/s]



[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
% Initialise matrices (N x N (x 9))
u = zeros(NXMAX,1);       % Velocity field
rho = ones(NXMAX,1);       % Density field
k=2*pi/50;x=1:50;
% rho(x)=1+0.01*(1-cos(k*x));
% u(x)=-0.01/sqrt(3)*(1-cos(k*x));
% rho(1:NXMAX/2)=3;rho(NXMAX/2+1:end)=1;

[feq] = Equilibrium(rho,scheme(:,2).',scheme(:,1).',u);
f=feq;%接下来是为了迭代初值

omega=3;


                          

figure
%% LBM loop
lambda=10;
W=2*pi/lambda;
omega=1;
tic
for t = 1 : T
    plot(rho);     
    pause(0.01)
    % Compute collision
    [f] = Collide(rho,scheme(:,2),scheme(:,1),u,cs,feq,NXMAX,f,omega,T0);
    % Stream populations
    f = Stream(f,NXMAX,scheme(:,1));
    % Find macroscopic quantities
    [rho,u] = Macroscopic(f,scheme(:,1));
    % Call visualisation
    rho(NXMAX/2)=1+0.01*(cos(k*t));
    u(:)=0.0; %模拟均匀流条件；%需要小于0.8

end
toc


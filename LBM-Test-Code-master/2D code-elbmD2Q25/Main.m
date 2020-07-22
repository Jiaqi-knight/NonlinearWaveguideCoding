%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 2D lattice scheme

%% Initialisation
clc
% Assumes D2Q25 and 25th trajectory is rest particle
clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
c_scale = [0 1];        % Velocity scale for plotting

% Time
T = 10000;        % Number of time steps

% Lattice dimensions
N = 1000;  % Number of x lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:N;

% Discrete velocities for D1Q5

[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(1,5);


% Initialise matrices (N x N (x 9))
f = zeros(N,5);       % Distribution function
feq = zeros(N,5);     % Equilibrium function
u = zeros(N,1);       % Velocity field
rho = zeros(N,1);       % Density field

% Forcing
forces_xy = zeros(N,N,2);   % Force field (Cartesian)
forces_i = zeros(N,N,9);    % Force field (Lattice)
F_extern = 0.00;           % External force (could be gravity)
F_dir = 'x';                % String indicating the direction of the 
                            % external force (either 'x' or 'y')

% Create typing matrix 'f' for fluid site
LatTyp = cell(N,N);
LatTyp(:,:) = {'f'};
% Boundary labels 'b'
% LatTyp(1,:) = {'b'}; LatTyp(N,:) = {'b'};
% LatTyp(2:N-1,1) = {'b'}; LatTyp(2:N-1,N) = {'b'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = 0;         % Initial Velocity Magnitude (x,y) (t = 0)
uref = .02;         % Reference velocity for Reynolds number
u(:) = u0;   % x-velocity field
tau=0.;

% Density
rho1 = 3.0;           % Initial Density (arbitrary)
rho2 = 1.0;           % Initial Density (arbitrary)
rho(1:N/2) = rho1;    % Uniform Field (arbitrary)
rho(N/2+1:end) = rho2;    % Uniform Field (arbitrary)

% Reynolds number
Re =   300;
lref = N;   % Reference length (in LUs) for Reynold number

% Other quantities
nu = uref * lref / Re;      % Kinematic viscosity from Reynolds
tau=(nu/cs^2)
omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
                                
% Information for user
disp(['Relaxation Time = ' num2str(1/omega)])
disp(['Kinematic Viscosity = ' num2str(nu)])
disp(['Reynolds Number = ' num2str(Re)])
disp('Running...')
figure
%% LBM loop
tic
for t = 1 : T
    plot(rho);

    pause(0.1)
    
% rho(1)=rho1;
% rho(end)=rho2;

    % Compute collision
    [f, feq] = Collide(rho,scheme(:,2),scheme(:,1),u,cs,feq,N,f,omega,T0);


    % Stream populations
    f = Stream(f,N,scheme(:,1));

    % Find macroscopic quantities
    [rho,u] = Macroscopic(f,scheme(:,1));
    
    % Call visualisation

end
toc
%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 2D lattice scheme

%% Initialisation

% Assumes D2Q249 and 49th trajectory is rest particle
clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
c_scale = [0.7 1.3];        % Velocity scale for plotting

% Time
T = 1000;        % Number of time steps
D=2;Q=7;
% Lattice dimensions
N = 100;  % Number of x lattice sites
M = 100;  % Number of y lattice sites

% Lattice site coordinates dx = 1
xl = 1:1:N; yl = 1:1:M;
[x,y] = meshgrid(xl,yl); % Grid of coordinates for plotting
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q^2);
                 

% Initialise matrices (N x M (x 9))
f = ones(M,N,Q^2)./Q^2;       % Distribution function
feq = zeros(M,N,Q^2);     % Equilibrium function
u = zeros(M,N,D);       % Velocity field
rho = zeros(M,N);       % Density field

                     % external force (either 'x' or 'y')

% Create typing matrix 'f' for fluid site
LatTyp = cell(M,N);
LatTyp(:,:) = {'f'};
% Boundary labels 'b'
% LatTyp(1,:) = {'b'}; LatTyp(M,:) = {'b'};
% LatTyp(2:M-1,1) = {'b'}; LatTyp(2:M-1,N) = {'b'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = [0 0];         % Initial Velocity Magnitude (x,y) (t = 0)
uref = .02;         % Reference velocity for Reynolds number
u(:,:,1) = u0(1);   % x-velocity field
u(:,:,2) = u0(2);   % y-velocity field

% Density
rho0 = 1;           % Initial Density (arbitrary)
rho(:,:) = rho0;    % Uniform Field (arbitrary)
% rho(50,50) = 0.01;    % Uniform Field (arbitrary)

% Reynolds number
Re = 100;
lref = M;   % Reference length (in LUs) for Reynold number

% Other quantities
cs = 1/sqrt(3);   % Sound speed on D2Q9 lattice
nu = uref * lref / Re;      % Kinematic viscosity from Reynolds
omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
  tau=nu/cs^2;                              % for unit time step
                                
% Information for user
disp(['Relaxation Time = ' num2str(tau)])
disp(['Kinematic Viscosity = ' num2str(nu)])
disp(['Reynolds Number = ' num2str(Re)])
disp('Running...')


lambda=30;
W=2*pi/lambda;

%% LBM loop
tic
for t = 1 : T
    % Compute collision
    [f, feq] = Collide(T0,rho,scheme(:,3),scheme(:,1:2),u,cs,feq,N,M,f,omega);

%     % Apply boundary conditions
%     f = Boundary(f,u,u0,rho,N,M);
    
    % Stream populations
    f = Stream(f,N,M,scheme(:,1:2));

    % Find macroscopic quantities
    [rho,u] = Macroscopic(f,scheme(:,1:2));
    
    % Call visualisation
    Visuals(x,y,u,rho,t,c_scale);
    rho(ceil(N/2),ceil(M/2)) = 2*sin(W*t);   

end
toc
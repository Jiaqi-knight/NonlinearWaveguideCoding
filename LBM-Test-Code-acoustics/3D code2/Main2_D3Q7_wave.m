%% Lattice Boltzmann LBGK Test Code
% Adrian Harwood, The University of Manchester, UK
% Last Updated: March 2016

%% Initialisation

% Assumes D3Q7 and 7th trajectory is rest particle
clear
set(0,'DefaultFigureWindowStyle','docked')
close all

% Visuals
out_every = 10;

% Time
T = 100;        % Number of time steps

% Lattice dimensions
N = 100;     % Number of x lattice sites
M = 2;     % Number of y lattice sites
K = 100;     % Number of z lattice sites

% Reynolds number
Re = 5;
lref = M;   % Reference length (in LUs) for Reynold number

% Lattice site coordinates dx = 1
zl = 1:1:N; yl = 1:1:M; xl = 1:1:K;
[x,y,z] = meshgrid(yl,xl,zl);   % Grid of coordinates for plotting

% Discrete velocities for D2Q9
xi = [   1,	-1,  0,  0,  0,  0,   0;...
        0,   0,  1, -1,  0,  0,   0;...
        0,   0,  0,  0,  1, -1,   0   ];
    
% Weights for D3Q7
w = [   1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/4   ];

% Initialise matrices ((7 x) K x M x N)
feq = zeros(7,K,M,N);      % Equilibrium function
J = zeros(3,K,M,N);         % Velocity field
% rho = zeros(K,M,N);         % Density field
P = zeros(K,M,N);         % Pressure field

% Forcing
forces_xy = zeros(3,K,M,N);     % Force field (Cartesian)
% forces_i = zeros(7,K,M,N);     % Force field (Lattice)
F_extern = 0.0001;           % External force (could be gravity)
F_dir = 'x';                % String indicating the direction of the 
                            % external force ('x' or 'y' or 'z')

% Create typing matrix 'f' for fluid site 'b' for boundary
LatTyp = cell(K,M,N);
LatTyp(:,:,:) = {'f'};

% Initialise macroscopic quanitities
% Velocity Field
u0 = [0 .03 0];       % Initial Velocity Magnitude (x,y,z) (t = 0)
uref = .03;         % Reference velocity for Reynolds number
J(1,:,:,:) = u0(1);   % x-velocity field
J(2,:,:,:) = u0(2);   % y-velocity field
J(3,:,:,:) = u0(3);   % z-velocity field

% Density
rho0 = 0;               % Initial Density (arbitrary)
P(:,:,:) = rho0;      % Uniform Field (arbitrary)


lambda=10;
W=2*pi/lambda;
A=5;



% Other quantities
c=0.5; %wave speed
cs = 1/2;             % Sound speed on D3Q7 lattice
% nu = uref * lref / Re;      % Kinematic viscosity from Reynolds
% omega = 1 / ( (nu/cs^2) + .5 ); % Relaxation frequency related to viscosity
%                                 % for unit time step
tau=1/2;
omega=1/tau;
% Initialise mesoscopic quantities
feq = Equilibrium(P,w,xi,J,c,cs,feq,N,M,K,0); %wjq-fix
f = feq;
% f(7,ceil(N/2),ceil(M/2),ceil(K/2)) = 0.01;                              
% % Information for user
% disp(['Relaxation Time = ' num2str(1/omega)])
% disp(['Kinematic Viscosity = ' num2str(nu)])
% disp(['Reynolds Number = ' num2str(Re)])
% disp('Running...')

% Visualise initial state
tsignal=[];
tsignal=Visuals(tsignal,N,M,K,x,y,z,J,P,0);


%% LBM loop
tic
for t = 1 : T
    
%     % Generate force vectors
%     [forces_i, forces_xy] = Force(forces_i,forces_xy,LatTyp,...
%         N,M,K,F_extern,F_dir,P,J,omega,w,cs,xi);
    P(ceil(N/2),ceil(M/2),ceil(K/2)) = A*sin(W*t);   
    % Compute collision
    [f, feq] = Collide(P,w,xi,J,c,cs,feq,N,M,K,f,omega);

    % Apply boundary conditions
%   f = Boundary(f,N,M,K,xi,LatTyp);
    
    % Stream populations
    f = Stream(f,N,M,K,xi);

    % Find macroscopic quantities
    [P,J] = Macroscopic(f,xi,LatTyp);
    
    % Call visualisation
%     if (mod(t,out_every) == 0)
        tsignal=Visuals(tsignal,N,M,K,x,y,z,J,P,t);
%     end
    
end
toc

     tsignal.cubes(solutime).v(1,:,:,:)=J(1,:,:,:);

%% theoretical
ex=reshape(tsignal.cubes(101).v(1,50,1,1:100),[1,100]);
th=A*besselj(0,W*2.5*[-49:50]);
figure 
plot(th./max(th));
hold on
plot(ex./max(ex))


 %% 整合数据，生成文件
%title=''; 
NAME = [date,'Duct'];  %存储文件夹   
NAME1 = ['Lattice_test','-',date];  %存储文件夹   
 
output_file_name=[NAME,'1.plt']; 
tsignal.Nvar=7;     
tsignal.varnames={'x','y','z','rho','ux','uy','uz'};
mat2tecplot(tsignal,output_file_name);
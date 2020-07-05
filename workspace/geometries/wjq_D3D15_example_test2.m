
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, August 2010
%%
%% The code does not take into account eny specific boundaru condiition.
%% Change to 3D by Jiaqi Wang，2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc
close all
subfunction_path1=genpath('C:\Users\Jiaqi-knight\Documents\GitHub\NonlinearWaveguideCoding\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\Jiaqi-knight\Documents\GitHub\NonlinearWaveguideCoding\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\Jiaqi-knight\Documents\GitHub\NonlinearWaveguideCoding\workspace\differential_geometry-master\matlab');
subfunction_path4=genpath('C:\Users\Jiaqi-knight\Documents\GitHub\NonlinearWaveguideCoding\workspace\Lattice_Boltzmann_method-master\matlab');
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
addpath(subfunction_path4);

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ny = 100;                    % Number of lines   (cells in the y direction)
Nx = 100;                    % Number of columns (cells in the x direction)
Nz = 100;                    % Number of columns (cells in the z direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
Lx = 0.5;                     % Maximun dimenssion in the x direction [m]
Ly = 0.5;                     % Maximun dimenssion on th y direction  [m]
Lz = 0.5;                     % Maximun dimenssion on th z direction  [m]
Dx = Lx/Nx;                    % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p;       % lattice time step


% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 1.985;                                    % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2
invcssq=1/cs2;
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity


% Block 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice properties for the D2Q9 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=3;Q=15; 
[scheme,cs,cssq,T0]=initializeELBMscheme(D,Q);
                                            % number of directions of the D3Q27 model
% C_x=scheme(:,1);                       % velocity vectors in x
% C_y=scheme(:,2);                       % velocity vectors in y
% C_z=scheme(:,3);                       % velocity vectors in y
% W  =scheme(:,4); 
C_x=[0  1  -1  0  0  0  0  1 -1  1 -1   1  -1  1 -1 ];
C_y=[0  0   0  1 -1  0  0  1 -1 -1  1  -1  1  -1  1 ];
C_z=[0  0   0  0  0  1 -1  1 -1 -1  1   1 -1  -1  1 ];
w0=2/9;w1=1/9;w3=1/72;
W  =[w0 w1 w1  w1 w1 w1 w1 w3 w3 w3 w3 w3 w3 w3 w3 ];
scheme(:,1)=C_x;
scheme(:,2)=C_y;
scheme(:,3)=C_z;
scheme(:,4)=W;


f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.


% Array of distribution and relaxation functions
f=  zeros(Nx,Ny,Nz,Q);
feq=zeros(Nx,Ny,Nz,Q);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/Q;
ux = zeros(Nx,Ny,Nz);
uy = zeros(Nx,Ny,Nz);
uz = zeros(Nx,Ny,Nz);

rho_l = 0.01;   % initial disturbance
f(Nx/2, Ny/2, Nz/2, 1) = rho_l;


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ta = 1 : 250*sqrt(3)
    
    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %硬壁面边界条件
f(:,       :,       :,      1)   =  f(  :,      :,       :,       1); %  0, 0, 0;
f(3:Nx,    :,       :,      2)   =  f( 2:Nx-1,  :,       :,       2); %  1, 0, 0;
f(1:Nx-2,  :,       :,      3)   =  f( 2:Nx-1,  :,       :,       3); % -1, 0, 0;
f(:,      3:Ny,     :,      4)   =  f(  :,     2:Ny-1,   :,       4); %  0, 1, 0;
f(:,      1:Ny-2,   :,      5)   =  f(  :,     2:Ny-1,   :,       5); %  0,-1, 0;
f(:,       :,      3:Nz,    6)   =  f(  :,      :,      2:Nz-1,   6); %  0, 0, 1;
f(:,       :,      1:Nz-2,  7)   =  f(  :,      :,      2:Nz-1,   7); %  0, 0,-1;
f(3:Nx,   3:Ny,    3:Nz,    8)   =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   8);%  1, 1, 1;
f(1:Nx-2, 1:Ny-2,  1:Nz-2,  9)   =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   9);% -1,-1,-1;
f(3:Nx,   3:Ny,    1:Nz-2,  10)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   10);%  1, 1,-1;
f(1:Nx-2, 1:Ny-2,  3:Nz,    11)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   11);% -1,-1, 1;
f(3:Nx,   1:Ny-2,  3:Nz,    12)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   12);%  1,-1, 1;
f(1:Nx-2, 3:Ny,    1:Nz-2,  13)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   13);% -1, 1,-1;
f(3:Nx,   1:Ny-2,  1:Nz-2,  14)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   14);%  1,-1,-1;
f(1:Nx-2, 3:Ny,    3:Nz,    15)  =  f( 2:Nx-1,  2:Ny-1, 2:Nz-1,   15);% -1, 1, 1;


% f(3:Nx,   3:Ny,     :,      8)   =  f( 2:Nx-1,  2:Ny-1,  :,       8); %  1, 1, 0;
% f(1:Nx-2, 1:Ny-2,   :,      9)   =  f( 2:Nx-1,  2:Ny-1,  :,       9); % -1,-1, 0;
% f(3:Nx,   1:Ny-2,   :,      10)  =  f( 2:Nx-1,  2:Ny-1,  :,       10);%  1,-1, 0;
% f(1:Nx-2, 3:Ny,     :,      11)  =  f( 2:Nx-1,  2:Ny-1,  :,       11);% -1, 1, 0;
% f(3:Nx,    :,      3:Nz,    12)  =  f( 2:Nx-1,  :,      2:Nz-1,   12);%  1, 0, 1;
% f(1:Nx-2,  :,      1:Nz-2,  13)  =  f( 2:Nx-1,   :,     2:Nz-1,   13);% -1, 0,-1;
% f(3:Nx,    :,      1:Nz-2,  14)  =  f( 2:Nx-1,  :,      2:Nz-1,   14);%  1, 0,-1;
% f(1:Nx-2,  :,      3:Nz,    15)  =  f( 2:Nx-1,   :,     2:Nz-1,   15);% -1, 0, 1;
% f(:,      3:Ny,    3:Nz,    16)  =  f(  :,     2:Ny-1,  2:Nz-1,   16);%  0, 1, 1;
% f(:,      1:Ny-2,  1:Nz-2,  17)  =  f(  :,     2:Ny-1,  2:Nz-1,   17);%  0,-1,-1;
% f(:,      3:Ny,    1:Nz-2,  18)  =  f(  :,     2:Ny-1,  2:Nz-1,   18);%  0, 1,-1;
% f(:,      1:Ny-2,  3:Nz,    19)  =  f(  :,     2:Ny-1,  2:Nz-1,   19);%  0,-1, 1;


    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,4);    
    

    
    % Determining the velocities according to Eq.() (see slides)
    
    for k=1:Q
    ux = ux+C_x(k).*f(:,:,k);
    uy = uy+C_y(k).*f(:,:,k);
    uz = uz+C_z(k).*f(:,:,k);
    end
    ux=ux./rho;
    uy=uy./rho;
    uz=uz./rho;
    f1=3.;
    f2=4.5;
    f3=1.5;  
    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uxsq=ux.^2;
    uysq=uy.^2;
    uzsq=uz.^2;
    usq=uxsq+uysq+uzsq;
%     rt0= w0*rho;
%     rt1= w1*rho;
%     rt2= w2*rho; 
for i=1:Q
    feq(:,:,:,i) = scheme(i,4).*rho.*(1 ...
    +(scheme(i,1).*ux+scheme(i,2).*uy+scheme(i,3).*uz).*f1 ...
    +((scheme(i,1).*ux+scheme(i,2).*uy+scheme(i,3).*uz).^2).*f2 ...
    -f3.*usq);   
end
    
    
    
%     feq(:,:,1)= rt1 .*(1 +f1*ux +f2.*uxsq -f3*usq);
%     feq(:,:,2)= rt1 .*(1 +f1*uy +f2*uysq -f3*usq);
%     feq(:,:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
%     feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
%     feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
%     feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
%     feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
%     feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
%     feq(:,:,9)= rt0 .*(1 - f3*usq);
%     
    
    
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f= (1-omega)*f + omega*feq;
    
    
    
    % Ploting the results in real time
    surf(rho(:,:,Nz/2)-1), view(2), shading flat, axis equal, caxis([-.00001 .00001])
    grid off
    pause(.0001)
    
end %  End main time Evolution Loop




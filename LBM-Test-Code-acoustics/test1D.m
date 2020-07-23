
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, August 2010
%%
%% The code does not take into account eny specific boundaru condiition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc
close all


% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 300;                    % Number of lines   (cells in the y direction)
c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
Lx = .5;                      % Maximun dimenssion in the x direction [m]
Dx = Lx/N                    % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p       % lattice time step


% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 1.9;                                      % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity


% Block 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice properties for the D2Q9 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c=3 ;                                             % number of directions of the D2Q9 model
C_x=[0 1  -1];                       % velocity vectors in x
w0=4/6. ; w1=1/6. ; w2=1/6.;                    % lattice weights
W = [w0 w1 w1];
f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.


% Array of distribution and relaxation functions
f=zeros(N,N_c);                                 
feq=zeros(N,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:)=rho_l/3;   
ux = zeros(N, 1);

rho_l = 0.01;   % initial disturbance
f(N/2,3) = rho_l;


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ta = 1 : 150*sqrt(3)
    
    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %这个貌似是反射边界条件
    f(:,2) = [f(1:2,2);f(2:N-1,2)];
    f(:,3) = [f(2:N-1,3);f(N-1:N,3)];
    
    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,2); 
    

        
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(2).*f(:,2)+C_x(3).*f(:,3))./rho ;
    mean(mean(mean(ux)))
    mean(mean(mean(f(:,:,:))))
    aola = 6        


    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    uxsq=ux.^2; 
    usq=uxsq; 
     
    feq(:,1)= rt0 .*(1 - f3*usq);
    feq(:,2)= rt1 .*(1 +f1*ux +f2*uxsq -f3*usq);
    feq(:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
%     feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
%     feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
%     feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
%     feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
%     feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
    
    
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f= (1-omega)*f + omega*feq; 
    
    
    
% Ploting the results in real time   
plot(rho-1), view(2)%, shading flat, axis equal, caxis([-.00001 .00001])
grid off
pause(.0001) 

end %  End main time Evolution Loop

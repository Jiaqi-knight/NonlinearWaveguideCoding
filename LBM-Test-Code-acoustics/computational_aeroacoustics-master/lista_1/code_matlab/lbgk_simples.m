
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

Nr = 300;                    % Number of lines   (cells in the y direction)
Mc = 300;                    % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
Lx = .5;                      % Maximun dimenssion in the x direction [m]
Ly = 0.0833;                  % Maximun dimenssion on th y direction  [m]
Dx = Lx/Mc                    % Lattice space (pitch)
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

N_c=9 ;                                             % number of directions of the D2Q9 model
C_x=[1 0 -1  0 1 -1 -1  1 0];                       % velocity vectors in x
C_y=[0 1  0 -1 1  1 -1 -1 0];                       % velocity vectors in y
w0=16/36. ; w1=4/36. ; w2=1/36.;                    % lattice weights
W = [w1 w1 w1 w1 w2 w2 w2 w2 w0];
f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.


% Array of distribution and relaxation functions
f=zeros(Nr,Mc,N_c);                                 
feq=zeros(Nr,Mc,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/9;   
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);

rho_l = 0.01;   % initial disturbance
f(Nr/2, Mc/2,9) = rho_l;


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ta = 1 : 150*sqrt(3)
    
    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
    f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
    f(:,:,3) = [f(:,2:Mc-1,3) f(:,Mc-1:Mc,3)];
    f(:,:,4) = [f(2:Nr-1,:,4);f(Nr-1:Nr,:,4)];
    f(:,:,5) = [f(:,1:2,5) f(:,2:Mc-1,5)];
    f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];
    f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)];
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)];
    
    

    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3); 
    

        
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
        

    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    uxsq=ux.^2; 
    uysq=uy.^2; 
    usq=uxsq+uysq; 
       
    feq(:,:,1)= rt1 .*(1 +f1*ux +f2.*uxsq -f3*usq);
    feq(:,:,2)= rt1 .*(1 +f1*uy +f2*uysq -f3*usq);
    feq(:,:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
    feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
    feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
    feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
    feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
    feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
    feq(:,:,9)= rt0 .*(1 - f3*usq);
    
    
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f= (1-omega)*f + omega*feq; 
    
    
    
% Ploting the results in real time   
surf(rho-1), view(2), shading flat, axis equal, caxis([-.00001 .00001])
grid off
pause(.0001) 

end %  End main time Evolution Loop




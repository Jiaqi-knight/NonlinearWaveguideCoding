%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flow past a square cylinder
%% By Andre M. N. Spillere, November 2015
%% Based on Andrey R. da Silva, August 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc
close all


plot_density = 0; % Plot animation? 0 = no / 1 = yes
plot_velocity = 0;

M = 0.10;     % Mach number
Nsquare = 1; % number of cells away from the body for FW-H surface

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = 500;                    % Number of lines   (cells in the y direction)
Mc = 500;                    % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
Lx = .5;                      % Maximun dimenssion in the x direction [m]
Dx = Lx/Mc;                    % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p;       % lattice time step

% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 1.93;                                      % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity
zeta = c_p/cs;

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



% Block 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ABC termination
a = 0.3;                    % sigma_m coefficient
D = 30;                     % number of cells in the anechoic region

% Full domain
Mc = Mc + 2*D;
Nr = Nr + 2*D;

stabilization_time = 0;
measurement_time = 2e5;
total_time = stabilization_time + measurement_time;


sigma_anechoic = zeros(Nr,Mc,N_c);

sigma = a*linspace(1/D,1,D).^2;
sigma_right = sigma;
sigma_left = fliplr(sigma);
sigma_up = sigma';
sigma_bottom = fliplr(sigma)';

sigma_anechoic(31:530,1:30,:) = repmat(sigma_left,500,1,N_c);
sigma_anechoic(31:530,531:560,:) = repmat(sigma_right,500,1,N_c);
sigma_anechoic(1:30,31:530,:) = repmat(sigma_bottom,1,500,N_c);
sigma_anechoic(531:560,31:530,:) = repmat(sigma_up,1,500,N_c);
sigma_anechoic(531:560,1:30,:) = (repmat(sigma_up,1,30,N_c).^2 + repmat(sigma_left,30,1,N_c).^2).^0.5;
sigma_anechoic(531:560,531:560,:) = (repmat(sigma_up,1,30,N_c).^2 + repmat(sigma_right,30,1,N_c).^2).^0.5;
sigma_anechoic(1:30,1:30,:) = (repmat(sigma_bottom,1,30,N_c).^2 + repmat(sigma_left,30,1,N_c).^2).^0.5;
sigma_anechoic(1:30,531:560,:) = (repmat(sigma_bottom,1,30,N_c).^2 + repmat(sigma_right,30,1,N_c).^2).^0.5;

% % Anechoic termination
ux_anechoic = M*cs + zeros(total_time,1);
uy_anechoic = zeros(total_time,1);
rho_anechoic = ones(total_time,1);
    
ux_anechoic_sq = ux_anechoic.^2;
uy_anechoic_sq = uy_anechoic.^2;
u_anechoic_sq = ux_anechoic_sq + uy_anechoic_sq;

% Wall vector
xl = [260 300 300 260 260 ];
yl = [300 300 260 260 300 ];

[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,xl,yl);

% Measurement points
x0 = xl(1) - Nsquare;
xf = xl(2) + Nsquare;
y0 = yl(3) - Nsquare;
yf = yl(2) + Nsquare;
vecx = x0:xf;
Nx = length(vecx);
vecy = y0:yf;
Ny = length(vecy);
square = [vecx(2:end), repelem(xf,Ny-1), fliplr(vecx(2:end)), repelem(x0,Ny-1); repelem(yf,Nx-1), fliplr(vecy(2:end)), repelem(y0,Nx-1), vecy(2:end)];
x_direction = square(1,:);
y_direction = square(2,:);

%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of distribution and relaxation functions
f=zeros(Nr,Mc,N_c);
feq=zeros(Nr,Mc,N_c);

% Target functions
f_target_anechoic = zeros(Nr,Mc,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/9;
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);

% Mean pressure and velocity at x1
pressure = zeros(1,total_time);
velocity = zeros(1,total_time);

% Time history
hist_rho = zeros(length(x_direction),total_time);
hist_ux = zeros(length(x_direction),total_time);
hist_uy = zeros(length(x_direction),total_time);
hist_u = zeros(length(x_direction),total_time);
hist_p = zeros(length(x_direction),total_time);


count = 1;
for ta = 1 : total_time
    
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
    
    
    %%%

    G = f;
    G(vec1) = f(vec3);
    G(vec2) = f(vec4);
    G(vec3) = f(vec1);
    G(vec4) = f(vec2);
    G(vec5) = f(vec7);
    G(vec6) = f(vec8);
    G(vec7) = f(vec5);
    G(vec8) = f(vec6);
    f = G;
   
    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3);
    p = rho*cs2;

    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;
     
    % Determining the velocities according to Eq.() (see slides)
    ux = 0;
    uy = 0;
    for i = 1:8
        ux = ux + C_x(i).*f(:,:,i);
        uy = uy + C_y(i).*f(:,:,i);
    end
    ux = ux./rho;
    uy = uy./rho;
    
    u = (ux.^2 + uy.^2).^0.5; 
    
%   Save pressure and velocity history
    for i = 1:length(x_direction)
        hist_rho(i,ta) = rho(y_direction(i),x_direction(i));
        hist_ux(i,ta) = ux(y_direction(i),x_direction(i));
        hist_uy(i,ta) = uy(y_direction(i),x_direction(i));
        hist_u(i,ta) = u(y_direction(i),x_direction(i));
        hist_p(i,ta) = p(y_direction(i),x_direction(i));
    end
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
    
    f_target_anechoic(:,:,1) = w1*rho_anechoic(ta).*(1 +f1*ux_anechoic(ta) +f2*ux_anechoic_sq(ta) -f3*u_anechoic_sq(ta));
    f_target_anechoic(:,:,2) = w1*rho_anechoic(ta).*(1 +f1*uy_anechoic(ta) +f2*uy_anechoic_sq(ta) -f3*u_anechoic_sq(ta));
    f_target_anechoic(:,:,3) = w1*rho_anechoic(ta).*(1 -f1*ux_anechoic(ta) +f2*ux_anechoic_sq(ta) -f3*u_anechoic_sq(ta));
    f_target_anechoic(:,:,4) = w1*rho_anechoic(ta).*(1 -f1*uy_anechoic(ta) +f2*uy_anechoic_sq(ta) -f3*u_anechoic_sq(ta));
    f_target_anechoic(:,:,5) = w2*rho_anechoic(ta).*(1 +f1*(+ux_anechoic(ta)+uy_anechoic(ta)) +f2*(+ux_anechoic(ta)+uy_anechoic(ta)).^2 -f3.*u_anechoic_sq(ta));
    f_target_anechoic(:,:,6) = w2*rho_anechoic(ta).*(1 +f1*(-ux_anechoic(ta)+uy_anechoic(ta)) +f2*(-ux_anechoic(ta)+uy_anechoic(ta)).^2 -f3.*u_anechoic_sq(ta));
    f_target_anechoic(:,:,7) = w2*rho_anechoic(ta).*(1 +f1*(-ux_anechoic(ta)-uy_anechoic(ta)) +f2*(-ux_anechoic(ta)-uy_anechoic(ta)).^2 -f3.*u_anechoic_sq(ta));
    f_target_anechoic(:,:,8) = w2*rho_anechoic(ta).*(1 +f1*(+ux_anechoic(ta)-uy_anechoic(ta)) +f2*(+ux_anechoic(ta)-uy_anechoic(ta)).^2 -f3.*u_anechoic_sq(ta));
    f_target_anechoic(:,:,9) = w0*rho_anechoic(ta).*(1 - f3*u_anechoic_sq(ta));    
      
    
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation)   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f = (1-omega)*f + omega*feq - sigma_anechoic.*(feq - f_target_anechoic);
    
    if rem(ta,1000) < 1
        disp(['Time step ' num2str(ta)])
        figure(1)
        imagesc(u) %caxis([0 M*cs])
        set(gca,'YDir','normal')
        line(xl,yl,'color','k')
        grid off
        title(['N = ' num2str(ta) ])
        pause(.001)
    end 
    
    %Ploting the results in real time
    if plot_density == true    
        surf(rho-1), view(2), shading flat, axis equal, caxis([-.001 .001])
        line(xl,yl,'color','k')
        grid off
        pause(.00001)
    end

    if plot_velocity == true
        imagesc(u) %caxis([0 M*cs])
        line(x_direction,y_direction,'color','r')    
        line(xl,yl,'color','k')   
        set(gca,'YDir','normal')
        grid off
        pause(.00001)
    end
end %  End main time Evolution Loop

save 'mach010_big.mat'
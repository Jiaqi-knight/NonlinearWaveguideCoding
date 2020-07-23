%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, August 2010
%%
%% The code does not take into account eny specific boundaru condiition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc
close all

tic
% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = 502;                    % Number of lines   (cells in the y direction)
Mc = 502;                    % Number of columns (cells in the x direction)

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
omega = 1.93;                                      % Relaxation frequency
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

% Calculando a condiÃ§Ã£o anecÃ³ica
% 4.0.1 - Adding conditions anechoic
distance = 30;
% condicao anecoica para cima
growth_delta = 0.5;
[sigma_mat9_cima Ft_cima] = build_anechoic_condition(Mc, ... 
Nr, distance, growth_delta);
% condicao anecoica para baixo
growth_delta = -0.5;
[sigma_mat9_baixo Ft_baixo] = build_anechoic_condition(Mc, ... 
Nr, distance, growth_delta);
% condicao anecoica para esquerda
growth_delta = -1;
[sigma_mat9_esquerda Ft_esquerda] = build_anechoic_condition(Mc, ... 
Nr, distance, growth_delta);
% condicao anecoica para direita
 growth_delta = 1;
 [sigma_mat9_direito Ft_direito] = build_anechoic_condition(Mc, ... 
 Nr, distance, growth_delta);

% Vendo pontos de barreira da xirizuda
xl = [2 2 501 501 2];
yl = [2 501 501 2 2];
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,xl,yl);
%xl = [250 250];
%yl = [1 50];
xl = [231 231 271 271 231];
yl = [231 271 271 231 231];
[vec1_duto,vec2_duto,vec3_duto,vec4_duto,vec5_duto,vec6_duto,vec7_duto,vec8_duto] = ... 
crossing3(Nr,Mc,xl,yl);

%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construindo chirp
total_time = 40*Mc*sqrt(3);
% times = 0 : total_time - 1;
% initial_frequency = 0;
% frequency_max_lattice = 1.8/(2*pi*20*sqrt(3));
% source_chirp = chirp(times, ... 
% initial_frequency, times(end), frequency_max_lattice);
% % vetores de pressao e velocidade de particula
struct_pressure = {};
Angle = (0:3:180)*pi/180;
raio = 210;
   for ponto = 1 : length(Angle)
       pressure(1:total_time) = 0;
       struct_pressure{ponto} = pressure(1:total_time);
   end

particle_velocity(1:total_time) = 0;
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

    G=f;
    f(vec1)=G(vec3);
    f(vec3)=G(vec1);
    f(vec2)=G(vec4);
    f(vec4)=G(vec2);
    f(vec5)=G(vec7);
    f(vec7)=G(vec5);
    f(vec6)=G(vec8);
    f(vec8)=G(vec6);

    G=f;
    f(vec1_duto)=G(vec3_duto);
    f(vec3_duto)=G(vec1_duto);
    f(vec2_duto)=G(vec4_duto);
    f(vec4_duto)=G(vec2_duto);
    f(vec5_duto)=G(vec7_duto);
    f(vec7_duto)=G(vec5_duto);
    f(vec6_duto)=G(vec8_duto);
    f(vec8_duto)=G(vec6_duto);

    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculando uma fonte ABC dentro do dominio
    %direita = 0.5
    Ma=[0.03 0.07 0.1];
    density_source = rho_l;% + 0.0001*source_chirp(ta);
    % Inlet
%     Ux_t = Ma(1)*cs;
%     Uy_t = 0;
%     point_y = 31;
%     distance_y = 469-31;
%     point_x = 2;
%     distance_x = 31;
%     direction = 0.5;
%     [sigma_source_inlet Ft_source_inlet] = build_source_anechoic(Nr, Mc, ...
%     density_source, Ux_t, Uy_t, point_y, ...
%     point_x, distance_x, distance_y, direction);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if ta == 1
    %rho(250, 150) = rho_l + 0.01*sin(ta);
    %end
    
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
    
      % colidindo tudo
    f = (1-omega)*f + omega*feq ...
    - sigma_mat9_cima.*(feq- Ft_cima) ...
    - sigma_mat9_direito.*(feq- Ft_direito)...
    - sigma_mat9_esquerda.*(feq- Ft_esquerda)...
    - sigma_mat9_baixo.*(feq- Ft_baixo);
    
    % Ploting the results in real time   
    
     %%
    if mod(ta, 30) == 0
        %vorticidade = curl(ux, uy);
        velocity = sqrt(ux.^2 + uy.^2);
        %imagesc(flip(vorticidade))
        disp('Progresso: ');
        disp((ta/total_time*100));
        imagesc(flip(velocity), [.00001 Ma(3)*cs])
        %axis equal
        pause(.000001)
    end
    %if ta>=total_time*0.8
    pressure_1(ta) = (rho(350, 250)-1)*cs2;
   % pressure_2(ta) = (rho(460, 200)-1)*cs2;
    %pressure_3(ta) = (rho(460, 300)-1)*cs2;
   % pressure_4(ta) = (rho(460, 400)-1)*cs2;
    vel_1(ta)=ux(350, 250);
   % vel_2(ta)=ux(460, 200);
   % vel_3(ta)=ux(460, 300);
   % vel_4(ta)=ux(460, 400);
   % end
   
   for ponto = 1 : length(Angle)
       pressure = struct_pressure{ponto};
       x = round(raio*cos(Angle(ponto)));
       y = round(raio*sin(Angle(ponto)));
       pressure(ta) = (rho(x, y)-1)*cs2;
       struct_pressure{ponto} = pressure;
   end
   
   
    
end %  End main time Evolution Loop
%%
figure(1)
 title('Pressao');
 imagesc(flip((rho-1)*cs), [-.001 .0022])
%%
 figure(2)
title('Magnitude de Velocidade');
 imagesc(flip(velocity),[.1 Ma(3)*cs])
toc
%%
 figure(3)
 vorticidade = curl(ux, uy);
title('Vorticidade');
 imagesc(flip(vorticidade),[-.0001 .001]);
toc
%%
figure(4)
 plot(pressure_1/pressure_1(end))
 hold on
 plot(vel_1/vel_1(end), 'r')
%%
figure(5)
p=-pressure_1(15000:end);
%plinha = p - abs(mean(p));
p_fft=abs(fft(p));
freq=linspace(0,1,length(p_fft));
st=40*freq/(Ma(3)*cs);
plot(st, p_fft)
axis([0.05 1.5 0 3])
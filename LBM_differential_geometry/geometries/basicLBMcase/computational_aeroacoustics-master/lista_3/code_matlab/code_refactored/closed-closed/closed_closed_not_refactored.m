
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
Lx = .198; 
Ly = 0.053;
Dx = 0.0005 

Nr = 53+4                    % Number of lines   (cells in the y direction)
Mc = 197                    % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
                    % Maximun dimenssion in the x direction [m]
                  % Maximun dimenssion on th y direction  [m]
                   % Lattice space (pitch)
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


 bcs=2;

  D_t=30;  % em número de celulas
  sigma_t=0.3;
  delta_t=0:D_t;




% Array of distribution and relaxation functions
f=zeros(Nr,Mc,N_c);                                 
feq=zeros(Nr,Mc,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/9;   
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);

%rho_l = 0.01;   % initial disturbance



%funcoes target - saida
Ux_t=0;
Uy_t=0;
U_t=Ux_t^2+Uy_t^2;
rho_t=rho_l;

coef1=  1/(2*cs2^2); %para uso na relaxacao
coef2= -1/(2*cs2);


Ft=zeros(Nr,Mc,9);
Fe=zeros(Nr,Mc,9);

Ft(:,:,9)= w0*rho_t.*(1+coef2*U_t);

Ft(:,:,1)= w1*rho_t.*(1 +Ux_t/cs2 +coef1*(Ux_t.^2 )+coef2*U_t);
Ft(:,:,2)= w1*rho_t.*(1 +Uy_t/cs2 +coef1*(Uy_t.^2 )+coef2*U_t);
Ft(:,:,3)= w1*rho_t.*(1 -Ux_t/cs2 +coef1*(Ux_t.^2 )+coef2*U_t);
Ft(:,:,4)= w1*rho_t.*(1 -Uy_t/cs2 +coef1*(Uy_t.^2 )+coef2*U_t);

Ft(:,:,5)= w2*rho_t.*(1 +(+Ux_t+Uy_t)/cs2 +coef1*((+Ux_t+Uy_t).^2) +coef2*U_t);
Ft(:,:,6)= w2*rho_t.*(1 +(-Ux_t+Uy_t)/cs2 +coef1*((-Ux_t+Uy_t).^2) +coef2*U_t);
Ft(:,:,7)= w2*rho_t.*(1 +(-Ux_t-Uy_t)/cs2 +coef1*((-Ux_t-Uy_t).^2) +coef2*U_t);
Ft(:,:,8)= w2*rho_t.*(1 +(+Ux_t-Uy_t)/cs2 +coef1*((+Ux_t-Uy_t).^2) +coef2*U_t);
% 




sigma=sigma_t*(delta_t/D_t).^2;
sigma_mat=[];
for i=1:Nr  % ver se tem jeito melhor de concatenar as matrizes
    sigma_mat=cat(1,sigma,sigma_mat);
end

sigmat=sigma_mat;

sigma_mat=[zeros(Nr,Mc-D_t-1) sigmat];
sigma_mat2=[sigmat zeros(Nr,Mc-D_t-1)];

% Condicao anecoica no final
inside= (round((53-6.3)/2)+3):(round((53+6.3)/2)+1);
outside=1:Nr;
for i=1:length(inside)  
    outside=outside(outside~=inside(i));
end
    
    
sigma_mat(outside ,:,:)=0;
sigma_mat2(outside ,:,:)=0;

sigma_mat9=[];
for i=1:9
sigma_mat9=cat(3,sigma_mat,sigma_mat9);
end

sigma_mat9e=[];
for i=1:9
sigma_mat9e=cat(3,sigma_mat2,sigma_mat9e);
end
% sigma_mat9e=zeros(Nr,Mc,9);
%    sigma_mat9e=fliplr(sigma_mat9);


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=0.001;
lambda=25;
freq=cs/lambda;
% REVER
xl=[   0           (30+61)     (30+61) (30+61+37) (30+61+37) (30+61+37+69)];
yl=[  (53-6.3)/2  (53-6.3)/2    0         0       (53-6.3)/2 (53-6.3)/2  ]+2;
% xl=[];
% yl=[];
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,xl,yl);
% 
yl2=[(53+6.3)/2 (53+6.3)/2    53      53 (53+6.3)/2  (53+6.3)/2  ]+2;
[vec12,vec22,vec32,vec42,vec52,vec62,vec72,vec82] = crossing3(Nr,Mc,xl,yl2);


taf= 2*Mc*sqrt(3);
figure;
for ta = 1 : taf
    
    
  
    
    
   
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
    f(vec12)=G(vec32);
    f(vec32)=G(vec12);
    f(vec22)=G(vec42);
    f(vec42)=G(vec22);
    f(vec52)=G(vec72);
    f(vec72)=G(vec52);
    f(vec62)=G(vec82);
    f(vec82)=G(vec62);
    

    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3); 
    
% %     Perturbacao
% %     rho(150,150)=rho_l+A*sin(2*pi*freq*(ta-1));
    if ta==1
    rho(inside,3)=1.005*rho_l;
    end 
    
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    % Determining dthe velocities according to Eq.() (see slides)    
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
    
    
    
    %funcoes target - saida


Ux_e=0;
Uy_e=0;
U_e=Ux_e^2+Uy_e^2;
lambda=25/ta;
freq=cs/lambda;
rho_e= rho_l+A*sin(2*pi*freq*(ta-1));

Fe(:,:,9)= w0*rho_e.*(1+coef2*U_e);

Fe(:,:,1)= w1*rho_e.*(1 +Ux_e/cs2 +coef1*(Ux_e.^2 )+coef2*U_e);
Fe(:,:,2)= w1*rho_e.*(1 +Uy_e/cs2 +coef1*(Uy_e.^2 )+coef2*U_e);
Fe(:,:,3)= w1*rho_e.*(1 -Ux_e/cs2 +coef1*(Ux_e.^2 )+coef2*U_e);
Fe(:,:,4)= w1*rho_e.*(1 -Uy_e/cs2 +coef1*(Uy_e.^2 )+coef2*U_e);

Fe(:,:,5)= w2*rho_e.*(1 +(+Ux_e+Uy_e)/cs2 +coef1*((+Ux_e+Uy_e).^2) +coef2*U_e);
Fe(:,:,6)= w2*rho_e.*(1 +(-Ux_e+Uy_e)/cs2 +coef1*((-Ux_e+Uy_e).^2) +coef2*U_e);
Fe(:,:,7)= w2*rho_e.*(1 +(-Ux_e-Uy_e)/cs2 +coef1*((-Ux_e-Uy_e).^2) +coef2*U_e);
Fe(:,:,8)= w2*rho_e.*(1 +(+Ux_e-Uy_e)/cs2 +coef1*((+Ux_e-Uy_e).^2) +coef2*U_e);
    
    
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%     f= (1-omega)*f + omega*feq; 
    
    f=omega*feq +(1-omega)*f-sigma_mat9.*(feq-Ft)-sigma_mat9e.*(feq-Fe) ;
    
% % Ploting the results in real time   
surf(rho-1), view(2), shading flat, axis equal,axis off, caxis([-.00001 .00001]), 
grid off
pause(.0000001) 
% % 

pre(:,ta)=cs*(rho(inside,32)-1);
prs(:,ta)=cs*(rho(inside,Mc-D_t-3)-1);

end %  End main time Evolution Loop

avg_prs=mean(prs);
avg_pre= mean(pre);

Prs=fft(avg_prs);
Pre=fft(avg_pre);

Z=Prs./Pre;
Znorm=Z*sqrt(3);

Zr=abs(real(Znorm));
Zi=abs(imag(Znorm));

f_util=2:size(Z,2)/2;
f_real=(c_p/cs)/Dx*f_util;

helm=f_real*0.072/c_p; % usar cs ou 340?

figure;
plot(f_real,Zr(2:end/2), 'r');
hold on
plot(f_real, Zi(2:end/2), 'b');

ylabel('Imped\E2ncia absoluta |Z_n_o_r_m|');
% xlabel('N\FAmero de Helmholtz');
xlabel('Frequencia [Hz]');
legend(' Real', 'Imagin\E1rio');
xlim([0 200e5]);


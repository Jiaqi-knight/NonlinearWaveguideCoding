%% Lattice Boltzmann LBGK Test Code
% Jiaqi Wang, SJTU, China
% Last Updated: March 2020
% vertification of the stability of 1D lattice scheme

%% Initialisation
clc
clear
%close all

viscosity=0.0001;%beta=1 %粘性
timeMax=200;
D=2;Q=9;
beta = 1/(2*viscosity+1); %eq-2.4
time = 1; %Start Time、
nx=300;
ny=300;
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
% Initialisation
populations=zeros(nx,ny,Q);
popequilibriums=zeros(nx,ny,Q);
populations(:,:,:)=1/Q;
ux = zeros(nx,ny);%Velocities at each lattice site are initialised as zero
uy = zeros(nx,ny); %Velocities at each lattice site are initialised as zero

k=2*pi/50;
populations(nx/2,ny/2,1)=2;




%function rho = LBMInitialise(nx)
rho = ones(nx,ny);
half = round(nx/2);
% rho(1:half) = 3;
% rho(half+1:end) = 1;
% rho(half,half)=0.01;

ux = zeros(nx,ny);%Velocities at each lattice site are initialised as zero
uy = zeros(nx,ny); %Velocities at each lattice site are initialised as zero




figure1 = figure;
axes1 = axes('Parent',figure1);
% hold(axes1,'on');
%% LBM loop


mov = moviein(timeMax);
while time <= timeMax
    %Stream
    populations = LBMPropagate(populations,nx,ny,scheme); %Propagate the populations
    %Macroscopic
    rho =sum(populations,3);
    ux = sum(repmat(permute(scheme(:,1).',[3,1,2]),nx,ny,1).*populations,3)./rho;
    uy = sum(repmat(permute(scheme(:,2).',[3,1,2]),nx,ny,1).*populations,3)./rho;
    
%     rho(nx/2,ny/2)=1+0.1*(cos(k*time));
    
    %Collide
    %case1:KBC
        popequilibriums(:,:,1) = rho*scheme(1,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy));
        popequilibriums(:,:,2) = rho*scheme(2,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux));
        popequilibriums(:,:,3) = rho*scheme(3,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1;
        popequilibriums(:,:,4) = rho*scheme(4,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
        popequilibriums(:,:,5) = rho*scheme(5,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
        popequilibriums(:,:,6) = rho*scheme(6,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
        popequilibriums(:,:,7) = rho*scheme(7,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
        popequilibriums(:,:,8) = rho*scheme(8,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
        popequilibriums(:,:,9) = rho*scheme(9,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
        [populations]=Kbc(beta,nx,ny,Q,rho,scheme,populations,popequilibriums,1);
    
    
    %case2:Qian
%     u2=ux.*ux+uy.*uy;
%     popequilibriums(:,:,1) = rho*scheme(1,3).*(1-u2*1.5);
%     popequilibriums(:,:,2) = rho*scheme(2,3).*(1-u2*1.5+3*ux+4.5*ux.^2);
%     popequilibriums(:,:,3) = rho*scheme(3,3).*(1-u2*1.5-3*ux+4.5*ux.^2);
%     popequilibriums(:,:,4) = rho*scheme(4,3).*(1-u2*1.5+3*uy+4.5*uy.^2);
%     popequilibriums(:,:,5) = rho*scheme(5,3).*(1-u2*1.5-3*uy+4.5*uy.^2);
%     popequilibriums(:,:,6) = rho*scheme(6,3).*(1-u2*1.5+3*(ux+uy)+4.5*(ux+uy).^2);
%     popequilibriums(:,:,7) = rho*scheme(7,3).*(1-u2*1.5+3*(ux-uy)+4.5*(ux-uy).^2);
%     popequilibriums(:,:,8) = rho*scheme(8,3).*(1-u2*1.5+3*(-ux+uy)+4.5*(-ux+uy).^2);
%     popequilibriums(:,:,9) = rho*scheme(9,3).*(1-u2*1.5+3*(-ux-uy)+4.5*(-ux-uy).^2);
%     populations_mirr=populations+2*(popequilibriums-populations);
%     populations = (1-beta)*populations +beta*populations_mirr;
    
    
    surf(rho-1), view(2), shading flat, axis equal%, caxis((-.01 .01)
    caxis([-0.01 0.01])
    colorbar(axes1);
    
    grid off% axis((0 nx -0.5 3.5));
    pause(0.0001)
    time = time + 1 %Increment time
end




















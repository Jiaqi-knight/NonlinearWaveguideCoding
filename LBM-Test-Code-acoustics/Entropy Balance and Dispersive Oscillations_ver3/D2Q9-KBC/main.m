clc
clear
close all
viscosity=0.4;%beta=1 %粘性
timeMax=100000;
D=2;Q=9;
beta = 1/(2*viscosity+1); %eq-2.4
time = 1; %Start Time、
nx=101;
ny=101;

%function rho = LBMInitialise(nx)
rho = ones(nx,ny);
half = round(nx/2);
% rho(1:half) = 3;
% rho(half+1:end) = 1;
% rho(half,half)=0.01;

ux = zeros(nx,ny);; %Velocities at each lattice site are initialised as zero
uy = zeros(nx,ny);; %Velocities at each lattice site are initialised as zero
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
%only for D1Q3--ELBM-平衡态
% populations(1,:) = rho./6.*( 3.*ux - 1 + 2*sqrt(1 + 3.*ux.^2));
% populations(2,:) = rho./6.*( -3.*ux - 1 + 2*sqrt(1 + 3.*ux.^2));
% populations(3,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*ux.^2));
% A_x=2- sqrt(1 + 3.*ux.^2);
% A_y=2- sqrt(1 + 3.*uy.^2);
% B_x=(2*ux+ sqrt(1 + 3.*ux.^2))./(1-ux);
% B_y=(2*uy+ sqrt(1 + 3.*uy.^2))./(1-uy);
% for k=1:Q
%     %Fabian-phd-34
%     populations(:,:,k) = rho.*scheme(k,3).*A_x.*A_y.*B_x.^scheme(k,1).*B_y.^scheme(k,2);
% end

populations(:,:,1) = rho*scheme(1,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy));
populations(:,:,2) = rho*scheme(2,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux));
populations(:,:,3) = rho*scheme(3,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1;
populations(:,:,4) = rho*scheme(4,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
populations(:,:,5) = rho*scheme(5,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
populations(:,:,6) = rho*scheme(6,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
populations(:,:,7) = rho*scheme(7,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
populations(:,:,8) = rho*scheme(8,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
populations(:,:,9) = rho*scheme(9,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;



% (populations,T) = LBMQuasiEquilibria(scheme,rho,ux,nx,T0); %Left moving,
% (populations,T)=entropyEquilibrium(nx,1,size(scheme,1),rho.',scheme(:,2).',scheme(:,1).',ux.',T0);%可能代码有bug,误差较大

figure1 = figure;
axes1 = axes('Parent',figure1);
% hold(axes1,'on');


mov = moviein(timeMax);
while time <= timeMax
    populations = LBMPropagate(populations,nx,ny,scheme); %Propagate the populations
    rho =sum(populations,3);
    ux = sum(repmat(permute(scheme(:,1).',[3,1,2]),nx,ny,1).*populations,3)./rho;
    uy = sum(repmat(permute(scheme(:,2).',[3,1,2]),nx,ny,1).*populations,3)./rho;

    
    % popequilibriums(1,:) = rho./6.*( 3.*ux - 1 + 2*sqrt(1 + 3.*ux.^2));
    % popequilibriums(2,:) = rho./6.*( -3.*ux - 1 + 2*sqrt(1 + 3.*ux.^2));
    % popequilibriums(3,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*ux.^2));
    %     A_x=2- sqrt(1 + 3.*ux.^2);
    %     A_y=2- sqrt(1 + 3.*uy.^2);
    %     B_x=(2*ux+ sqrt(1 + 3.*ux.^2))./(1-ux);
    %     B_y=(2*uy+ sqrt(1 + 3.*uy.^2))./(1-uy);
    %     for k=1:Q
    %         %Fabian-phd-34
    %         popequilibriums(:,:,k) = rho.*scheme(k,3).*A_x.*A_y.*B_x.^scheme(k,1).*B_y.^scheme(k,2);
    %     end
    %
    %     bot=1.6;top=2.4;err=1E-5;
    %     (alpha)=erfenfa(scheme,populations,popequilibriums,bot,top,err);
    
popequilibriums(:,:,1) = rho*scheme(1,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy));
popequilibriums(:,:,2) = rho*scheme(2,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux));
popequilibriums(:,:,3) = rho*scheme(3,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1;
popequilibriums(:,:,4) = rho*scheme(4,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
popequilibriums(:,:,5) = rho*scheme(5,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
popequilibriums(:,:,6) = rho*scheme(6,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
popequilibriums(:,:,7) = rho*scheme(7,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;
popequilibriums(:,:,8) = rho*scheme(8,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy));
popequilibriums(:,:,9) = rho*scheme(9,3).*(2-sqrt(1+3*ux.*ux)).*(2-sqrt(1+3*uy.*uy)).*((2*ux+sqrt(1+3*ux.*ux))./(1-ux)).^-1 .*((2*uy+sqrt(1+3*uy.*uy))./(1-uy)).^-1;

rho(half,half)=0.02*cos(2*pi/50*time);

[populations]=Kbc(beta,nx,ny,Q,rho,scheme,populations,popequilibriums,0);
    
    
    
    
    
    
    
    %     populations_mirr=populations+repmat(alpha,1,1,Q).*(popequilibriums-populations);
    %     populations = (1-beta)*populations +beta*populations_mirr;
    
    surf(rho-1), view(2), shading flat, axis equal%, caxis((-.01 .01)
%     caxis([-.001 .001])
    colorbar(axes1);
    
    grid off% axis((0 nx -0.5 3.5));
    pause(0.0001)
    time = time + 1 %Increment time
end




















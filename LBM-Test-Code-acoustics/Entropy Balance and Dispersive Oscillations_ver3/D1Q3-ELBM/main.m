clc
clear
close all
viscosity=0.1;%beta=1 %粘性
timeMax=100000;
D=1;Q=3;
beta = 1/(2*viscosity+1); %eq-2.4
time = 1; %Start Time、
nx=1000;
%function rho = LBMInitialise(nx)
%case1:
% rho = zeros(1,nx);
% half = floor(nx/2);
% rho(1:half) = 3;
% rho(half+1:end) = 1;
%case2:
rho = zeros(1,nx);
half = floor(nx/2);
rho(:) = 1;
rho(half) = 0.01;




u = 0; %Velocities at each lattice site are initialised as zero
[scheme,cs,cssq,invcs,invcssq,T0]=initializeELBMscheme(D,Q);
%only for D1Q3--ELBM-平衡态
populations(1,:) = rho./6.*( 3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
populations(2,:) = rho./6.*( -3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
populations(3,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*u.^2));
% [populations,T] = LBMQuasiEquilibria(scheme,rho,u,nx,T0); %Left moving,
% [populations,T]=entropyEquilibrium(nx,1,size(scheme,1),rho.',scheme(:,2).',scheme(:,1).',u.',T0);%可能代码有bug,误差较大



mov = moviein(timeMax);
x = 1:nx;
while time <= timeMax
populations = LBMPropagate(populations,nx,scheme); %Propagate the populations
rho =sum(populations,1);
u = scheme(:,1).'*populations./rho;

popequilibriums(1,:) = rho./6.*( 3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
popequilibriums(2,:) = rho./6.*( -3.*u - 1 + 2*sqrt(1 + 3.*u.^2));
popequilibriums(3,:) = 2.*rho./3.*(2 - sqrt(1 + 3.*u.^2));
% [popequilibriums,T] =entropyEquilibrium(nx,1,size(scheme,1),rho.',scheme(:,2).',scheme(:,1).',u.',T);


bot=1.6;top=2.4;err=1E-5;
[alpha]=erfenfa(scheme,populations,popequilibriums,bot,top,err);
populations_mirr=populations+repmat(alpha,Q,1).*(popequilibriums-populations);
populations = (1-beta)*populations +beta*populations_mirr;

plot(rho);
axis([0 nx -0.5 3.5]);
pause(0.0001)
time = time + 1; %Increment time
 
end




















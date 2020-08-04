%LBM Lattice Boltzmann Method
% LBM(MC,V,NX,TM,MOV,ESPIL,LIMIT,NORM) solves the shock tube problem
% with number of lattice
% sites NX with separation of one and viscosity V upto TM time steps of
% length one using method
% MC:
% 0 − LBGK polynomial equilibria,
% 1 − LBGK entropic equilibria,
% 2 − ELBM Newton iterations,
% 3 − ELBM Bisection method
% and entropic limiter
% LIMIT:
% 0 − No limiter
% 1 − Median filtering
% Output of a movie of the simulation can be controlled with
% MOV:
% 0 − No movie,
% 1 − With movie
% The accuracy of the root finding for ELBM can be controlled with
% EPSIL in all cases the root found will result in entropy production.
% Convergence is partly based on | | feq − f|| , the choice of the norm is
% controlled with
% NORM:
% 0 − L1 norm
% 1 − Entropic Norm;
function [rho, mov, convergence, relativeFiltering, alpha] = ...
LBM(methodChoice,viscosity,nx,timeMax,MOV,ELBMEpsilon,Limiter,Norm)
% Over−relaxation parameter
beta = 1/(2*viscosity+1); %eq-2.4
time = 1; %Start Time
rho = LBMInitialise(nx); %Densities at each lattice site are intialised as
%the standard shock tube problem
u = 0; %Velocities at each lattice site are initialised as zero
T=1/3; %初始温度
[populations,T] = LBMQuasiEquilibria(rho,u,methodChoice,nx,T); %Left moving,
%central and right moving populations at each lattice site respectively are
%initialised as the quasiequilibria
relativeFiltering = zeros(1,251);
if (MOV == 1)
mov = moviein(timeMax);
x = 1:nx;
else
mov = 0;
end
while time <= timeMax
populations = LBMPropagate(populations,nx); %Propagate the populations
% k=2*pi/50;
% populations(2,nx/2)=populations(2,nx/2)+0.1*(cos(k*time));
%keyboard
[rho, u] = LBMLatticeParameters(populations); %Density and velocity are
max(abs(u))
%calculated using the populations
[popequilibriums,T] = LBMQuasiEquilibria(rho,u,methodChoice,nx,T); %New
% quasiequilibria are found using new density and velocity values
[alpha , convergence] = LBMEntropicParameter(populations, ...
popequilibriums,methodChoice,nx,ELBMEpsilon,Norm);
% Finding the non trivial root for constant entropy in the ELBM, for
% normal LBGK alpha = 0
[limiterSites,shockRelative] = LBMLimiterSites(populations, ...
popequilibriums,nx,Limiter);
relativeFiltering(250 - shockRelative) = ...
relativeFiltering(250 - shockRelative) + 1;
populations = LBMCollide(populations,popequilibriums,beta,alpha,nx, ...
Limiter,limiterSites);
% Populations are collided with the quasiequilibria
if(~isreal(populations))
keyboard
return
end
if(MOV == 1) % If movie parameter is enabled record a frame
unlimitedSites = setdiff(1:nx,limiterSites);
plot(unlimitedSites,rho(unlimitedSites),'.',limiterSites, ...
rho(limiterSites),'r*');
axis([0 nx -0.5 3.5])
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            mov(:,time) = getframe;
end
time = time + 1; %Increment time
end


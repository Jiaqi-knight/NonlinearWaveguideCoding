%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the stationary analytical solution of the viscous wave
% equation for an outgoing cylindrical wave radiated from an infinite line
% source(see Eq. (2) from class exercise 1). 
%          
%                   p = cylin_wave(freq,nu,cs,A,x,phi)
% 
% 'freq' is the source frequency in Hertz;
% 'nu' is the physical kinematic viscosity of the medium [m^2/s];
% 'cs' is the physical speed of sound [m/s];
% 'A' is the amplitude of the pressure disturbance at the source [Pa];
% 'x' is the spatial coordinate associated with p [m];
% 'phi' is the phase angle of the source [rad]; 
% 
% The function returns the radial space history of the pressure p(x) along the
% spatial coordinate vector 'x'. 
%
% By Andrey R. da Silva, 2010
% Federal University of Santa Catarina
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p,x] = cylin_wave(freq,nu,cs,A,x,phi)



tau = 1/cs^2 * (2 * nu);
omega = 2*pi*freq;
alpha = omega/(cs*sqrt(2))*sqrt((sqrt(1+(omega*tau)^2)-1)/(1+(omega*tau)^2));
k = omega/cs - i*alpha;
x(1) = x(1)+eps;
p = real(A*exp(i*(phi-pi/2))*besselh(0,2,k*x));




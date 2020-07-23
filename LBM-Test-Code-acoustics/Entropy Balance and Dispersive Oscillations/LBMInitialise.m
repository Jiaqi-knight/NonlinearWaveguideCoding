%A function to intialise the lattice densities for the shock tube problem
function rho = LBMInitialise(nx)
rho = zeros(1,nx);
half = floor(nx/2);
rho(:) = 1;
%  rho(half+1) = 0.5;

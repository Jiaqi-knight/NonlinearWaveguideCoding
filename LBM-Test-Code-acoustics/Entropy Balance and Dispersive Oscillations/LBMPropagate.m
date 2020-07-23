%A function to propagate the populations
function populations = LBMPropagate(oldpopulations,nx)
populations = oldpopulations;
left = populations(1,1);
right = populations(3,nx);
populations(1,1:nx-1) = populations(1,2:nx);
populations(3,2:nx) = populations(3,1:nx-1);
populations(3,1) = left;
populations(1,nx) = right;

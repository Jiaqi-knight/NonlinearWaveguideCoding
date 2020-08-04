%A function to find density and velocity at each lattice site from the
%populations
function [rho,u] = LBMLatticeParameters(populations,scheme)
rho =sum(populations,1);
u = scheme(:,1).'*populations./rho;



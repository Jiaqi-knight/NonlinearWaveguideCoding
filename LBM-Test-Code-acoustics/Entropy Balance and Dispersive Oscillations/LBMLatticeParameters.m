%A function to find density and velocity at each lattice site from the
%populations
function [rho,u] = LBMLatticeParameters(populations)
rho = populations(1,:) + populations(2,:) + populations(3,:);
u = (populations(3,:) - populations(1,:))./rho;
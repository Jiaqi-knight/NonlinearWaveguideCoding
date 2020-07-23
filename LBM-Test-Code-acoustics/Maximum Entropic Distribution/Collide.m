function [f] = Collide(NXMAX,D,Q,T0,rho,w,c,u,cs,feq,N,f,omega)
% Loop through the lattice points to compute the new ditribution functions.
% Equilibrium based on:
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

% [feq] = Equilibrium(rho,w.',c.',u);
[feq] = entropyEquilibrium(NXMAX,D,Q,T0,rho,w.',c.',u);
% sum(feq,2)
% Recompute distribution function f
% fnew = (omega * (feq-f  ) ) + f ;
f=f-omega*(f-feq);%collision

end
    

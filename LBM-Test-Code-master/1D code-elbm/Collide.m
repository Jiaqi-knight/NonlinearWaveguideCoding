function [fnew, feq] = Collide(rho,w,c,u,cs,feq,N,f,omega,T0)
% Loop through the lattice points to compute the new ditribution functions.
% Equilibrium based on:
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

for i = 1:N  
        for k = 1:5
            % Compute f^eq
            feq(i,k) = rho(i) * w(k) * ...
                (1 + c(k)*u(i)/T0 +u(i)^2/2/T0^2*(c(k)^2-T0)+u(i)^3/6/T0^3*c(k)*(c(k)^2-3*T0));     
        end
        
end
% sum(feq,2)
% Recompute distribution function f
fnew = (omega * (feq-f  ) ) + f ;
    
end
    

function [ feq] = Equilibrium(rho,w,c,u)

feq=rho*w+3*(rho*w).*(u*c)+9/2*(rho*w).*(u*c).^2-3/2*(rho.*(u.^2))*w;

end


function [f] = Collide(rho,w,c,u,f,omega)
% Loop through the lattice points to compute the new ditribution functions.
% Equilibrium based on:
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

            feq = rho * w.' .* ...
                ( 1 + ...
                (3)*u*c.' + ...
                (9/2)*(u*c.').^2 - ...
                (3/2)* u.^2 );
%     rt0= w(1)*rho;
%     rt1= w(2)*rho;
%     rt2= w(3)*rho;
%     ux = u;
%     uxsq=ux.^2; 
%     usq=uxsq; 
%     
%     f1=3.;
%     f2=4.5;
%     f3=1.5;    
%     feq(:,1)= rt0 .*(1 - f3*usq);
%     feq(:,2)= rt1 .*(1 +f1*ux +f2*uxsq -f3*usq);
%     feq(:,3)= rt2 .*(1 -f1*ux +f2*uxsq -f3*usq);
% 

    f= (1-omega)*f + omega*feq; 

end
    

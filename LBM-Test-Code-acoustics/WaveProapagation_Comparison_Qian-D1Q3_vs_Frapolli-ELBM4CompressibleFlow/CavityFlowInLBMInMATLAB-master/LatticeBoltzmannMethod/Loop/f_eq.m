function feq  = f_eq(rho, u )
% This function is to compute the equilibrium distribution function of
% input is field varibles rho(1), u(2) Here rho and u is local grid node variable 
% Output is varibles feq(Q) equilibrium density function of the grid node
% needed global constant e(2,9), omega(1,9), Q, Nx, Ny, 
global Dimension Q ev omega c
feq = zeros(Q,1);
% meso variables
c_square = c^2. ;
a1 = 3. / c_square ; 
a2 = 9./2. / c_square^2. ;
a3 =  3./2. / c_square ;

edotu = zeros(Q,1);
for k = 1:Q
    for m = 1:Dimension
        edotu(k) = edotu(k) + ev(m,k)*u(m);
    end
end

u_square = 0;
for m = 1:Dimension
     u_square = u_square + u(m).^2; 
end

% compute f_eq
for k = 1:Q
   feq(k) = rho * omega(k) * (1 + a1.*edotu(k) + a2.*edotu(k).^2 - a3.*u_square );
end

end

% time aeq to O(Q*Dimension)


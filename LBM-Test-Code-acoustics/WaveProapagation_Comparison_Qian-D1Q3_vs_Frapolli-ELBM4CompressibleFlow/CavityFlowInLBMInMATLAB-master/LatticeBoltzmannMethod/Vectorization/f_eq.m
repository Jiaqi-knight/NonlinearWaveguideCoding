function [feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8,feq9]  = f_eq(rho, ux, uy)
% This function is to compute the equilibrium distribution function of
% whole flow field
% input is field varibles rho(1), ux(1), uy(1) Here rho and u is local grid node variable 
% Output is varibles feq(Q) equilibrium density function of the grid node
% needed global constant e(2,9), omega(1,9), Q, Nx, Ny, 
global  Q ex ey omega c
feq = zeros(Q,1);
% meso variables
c_square = c^2. ;
a1 = 3. / c_square ; 
a2 = 9./2. / c_square^2. ;
a3 =  3./2. / c_square ;

edotu = zeros(Q,1);
edotu(:) = ex(:)*ux + ey(:)*uy;          % O(Q)
u_square = ux^2 + uy^2; 

% compute f_eq
feq(:) = rho * omega(:) .* (1 + a1*edotu(:) + a2*edotu(:).^2 - a3*u_square );  % O(Q)

feq1=feq(1); feq2=feq(2); feq3=feq(3); feq4=feq(4); feq5=feq(5); 
feq6=feq(6); feq7=feq(7); feq8=feq(8); feq9=feq(9); 


end

% time aeq to O(Q*Dimension)


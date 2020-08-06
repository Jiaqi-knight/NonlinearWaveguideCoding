%% Overload of collision function to allow calculation of feq only for initialisation
function [feq]= LBM_collide( i, j, k, v, r ) 

% 	/* LBGK equilibrium function is represented as:
% 		feq_i = rho * w_i * ( 1 + u_a c_ia / cs^2 + Q_iab u_a u_b / 2*cs^4 )
% 	where
% 		Q_iab = c_ia c_ib - cs^2 * delta_ab
% 	and
% 		delta_ab is the Kronecker delta.
% 	*/
	
	%% Declare single feq value and intermediate values A and B
% 	 feq, A, B;
    global Grids dims c w cs

	%% Other declarations
	N_lim = length(Grids{r}.XPos);
	M_lim = length(Grids{r}.YPos);
	K_lim = length(Grids{r}.ZPos);
	
	%% Compute the parts of the expansion for feq

		%% 3D case -- IF CONTIGUOUS ONLY NEED 1 INDEX AND USE POINTER ARITHMETIC
		idx0 = idxmap(i,j,k,1,M_lim,K_lim,dims);
		idx1 = idxmap(i,j,k,2,M_lim,K_lim,dims);
		idx2 = idxmap(i,j,k,3,M_lim,K_lim,dims);

if (dims == 3)
		%% Compute c_ia * u_a which is actually the dot product of c and u
		A = (c(1,v) * Grids{r}.u(idx0)) + (c(2,v) * Grids{r}.u(idx1)) + (c(3,v) * Grids{r}.u(idx2));

% 		/*
% 		Compute second term in the expansion
% 		Q_iab u_a u_b = 
% 		(c_x^2 - cs^2)u_x^2 + (c_y^2 - cs^2)u_y^2 + (c_z^2 - cs^2)u_z^2
% 		+ 2c_x c_y u_x u_y + 2c_x c_z u_x u_z + 2c_y c_z u_y u_z
% 		*/

		B =	(power(c(1,v),2) - power(cs,2)) * power(Grids{r}.u(idx0),2) + ...
			(power(c(2,v),2) - power(cs,2)) * power(Grids{r}.u(idx1),2) + ...
			(power(c(3,v),2) - power(cs,2)) * power(Grids{r}.u(idx2),2) + ...
			2 * c(1,v)*c(2,v) * Grids{r}.u(idx0) * Grids{r}.u(idx1) + ...
			2 * c(1,v)*c(3,v) * Grids{r}.u(idx0) * Grids{r}.u(idx2) + ...
			2 * c(2,v)*c(3,v) * Grids{r}.u(idx1) * Grids{r}.u(idx2);
else
		%% 2D versions of the above
		A = (c(1,v) * Grids{r}.u(idx0)) + (c(2,v) * Grids{r}.u(idx1));

		B =	(power(c(1,v),2) - power(cs,2)) * power(Grids{r}.u(idx0),2) + ...
			(power(c(2,v),2) - power(cs,2)) * power(Grids{r}.u(idx1),2) + ...
			2 * c(1,v)*c(2,v) * Grids{r}.u(idx0) * Grids{r}.u(idx1);
end
	
	
	%% Compute f^eq
	idxijk = idxmap(i,j,k,M_lim,K_lim);
	feq = Grids{r}.rho(idxijk) * w(v) * (1 + (A / power(cs,2)) + (B / (2*power(cs,4))));



end

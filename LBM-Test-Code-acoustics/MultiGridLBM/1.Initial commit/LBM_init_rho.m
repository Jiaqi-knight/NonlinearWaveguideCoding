%% Initialise density on level r
function  LBM_init_rho (r) 
    global Grids dims cs rho_in u_0x u_0y u_0z a_x a_y a_z b_x b_y b_z kn km kk
	%% Get grid sizes
	 N_lim = length(Grids{r}.XPos);
	 M_lim = length(Grids{r}.YPos);
	 K_lim = length(Grids{r}.ZPos);

	for ( i = 1 : N_lim) 
		for ( j = 1 : M_lim) 	
			for ( k = 1 :K_lim) 

				%% Max velocity
				 u_in = [u_0x, u_0y, u_0z];

				%% Wave numbers
				 Lx = b_x - a_x;
				 Ly = b_y - a_y;
				 Lz = b_z - a_z;
				 k1 = 2*pi*kn / Lx;
				 k2 = 2*pi*km / Ly;
				 k3 = 2*pi*kk / Lz;

				 idx = idxmap(i,j,k,M_lim,K_lim);

				%% Only a 2D expression used here for now...
				Grids{r}.rho(idx) = rho_in - ...
					(1/power(cs,2)) * .25 * power(norm(u_in,2),2) * ...
					(cos(2*k1*Grids{r}.XPos(i))+cos(2*k2*Grids{r}.YPos(j)));
            end
        end
    end

%     end

% % Initialise velocity on level r
function LBM_init_vel(r)
    global Grids dims u_0x u_0y u_0z a_x a_y a_z b_x b_y b_z kn km kk
	% Max velocity
	u_in = [u_0x, u_0y, u_0z];
	
	% Wave numbers
	Lx = b_x - a_x;
	Ly = b_y - a_y;
	Lz = b_z - a_z;
	k1 = 2*pi*kn / Lx;
	k2 = 2*pi*km / Ly;
	k3 = 2*pi*kk / Lz;

	% Get grid sizes
	N_lim = length(Grids{r}.XPos);
	M_lim = length(Grids{r}.YPos);
	K_lim = length(Grids{r}.ZPos);

	for (i = 1 : N_lim) 
		for (j = 1 : M_lim) 
			for (k = 1 : K_lim) 


				% Only a 2D expression used here for now...
				idx = idxmap(i,j,k,1,M_lim,K_lim,dims);
				Grids{r}.u(idx) = -norm(u_in,2) * ...
					cos(k1*Grids{r}.XPos(i)) * sin(k2*Grids{r}.YPos(j));

				idx = idxmap(i,j,k,2,M_lim,K_lim,dims);
				Grids{r}.u(idx) = norm(u_in,2) * ...
					(k1/k2) * sin(k1*Grids{r}.XPos(i)) * cos(k2*Grids{r}.YPos(j));
            end
        end
    end
end

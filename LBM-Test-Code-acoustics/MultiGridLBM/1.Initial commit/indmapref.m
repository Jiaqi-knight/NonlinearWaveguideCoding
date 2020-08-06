% Routine to map the index of a coarse grid site to a corresponding fine site on the level below
function [fine_ind]=indmapref( coarse_i,  x_start,  coarse_j,  y_start,  coarse_k,  z_start) 

	% Initialise result
	
	% Map indices
	fine_ind(1)= 2*(coarse_i - x_start + 1) - 1 ;
	fine_ind(2)= 2*(coarse_j - y_start + 1) - 1 ;
	fine_ind(3)= 2*(coarse_k - z_start + 1) - 1 ;

end
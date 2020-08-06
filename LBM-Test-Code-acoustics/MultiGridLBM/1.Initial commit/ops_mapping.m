% /* File contains all the routines necessary to perform 
% mapping of positions or indices.
% */

% #include "stdafx.h"
% #include "LBM_globalvars.h"

% using namespace std;




     


% ***************************************************************************************************

% Routine to map the position of a coarse grid site to a corresponding fine site on the level below
function [fine_pos]=posmapref ( coarse_pos, fine_level,  direction,  plusminus) 

	% Mapping routine assuming refinement level of 2.
	% Returns the position in the level 0 reference frame of the first
	% element of the corresponding pair of 2 nodes in the adjacent finer grid.

	% Spacing
% 	 spacing;
	if (direction == 'x') 
		spacing = Grids{fine_level}.dx/2;
	 elseif (direction == 'y') 
		spacing = Grids{fine_level}.dy/2;
	 elseif (direction == 'z') 
		spacing = Grids{fine_level}.dz/2;
    end

	% Position
% 	 fine_pos;
	if (plusminus == '+') 
		fine_pos = coarse_pos + spacing;
	 elseif (plusminus == '-') 
		fine_pos = coarse_pos - spacing;
    end

    
end

% ***************************************************************************************************

% Routine to map the index of a coarse grid site to a corresponding fine site on the level below
function [fine_ind]=indmapref(coarse_i, x_start, coarse_j, y_start, coarse_k, z_start) 

	% Initialise result
% 	vector<int> fine_ind;
	
	% Map indices
	fine_ind.insert(fine_ind.begin(), 2*(coarse_i - x_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 1, 2*(coarse_j - y_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 2, 2*(coarse_k - z_start + 1) - 2 );

end
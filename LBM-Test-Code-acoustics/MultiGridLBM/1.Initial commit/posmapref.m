%% Routine to map the position of a coarse grid site to a corresponding fine site on the level below
function [fine_pos]=posmapref ( coarse_pos,  fine_level,  direction,  plusminus) 

	%% Mapping routine assuming refinement level of 2.
	%% Returns the position in the level 0 reference frame of the first
	%% element of the corresponding pair of 2 nodes in the adjacent finer grid.
     global Grids
	%% Spacing
% 	double spacing;
	if (direction == 'x') 
		spacing = Grids{fine_level}.dx/2;
    elseif (direction == 'y') 
		spacing = Grids{fine_level}.dy/2;
	elseif (direction == 'z') 
		spacing = Grids{fine_level}.dz/2;
    end

	%% Position
% 	double fine_pos;
	if (plusminus == '+') 
		fine_pos = coarse_pos + spacing;
    elseif (plusminus == '-') 
		fine_pos = coarse_pos - spacing;
    end

end

%% ***************************************************************************************************
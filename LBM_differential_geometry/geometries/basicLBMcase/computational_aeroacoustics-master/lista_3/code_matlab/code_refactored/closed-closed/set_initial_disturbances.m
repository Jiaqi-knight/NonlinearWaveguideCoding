function lattice = set_initial_disturbances(lattice, initial_disturbance_density, points_lattice)
	
	lattice_distribution = lattice{1};
	size_lattice = size(lattice_distribution(:, :, 1));
	number_lines_lattice = size_lattice(1);
	number_columns_lattice = size_lattice(2);

	horizontal_points = points_lattice{2};
	vertical_points = points_lattice{1};

	if length(horizontal_points) ~= length(vertical_points)
		disp('Quantity points not equal in X and Y axis.');
	else
		quantity_points = length(horizontal_points);
		for point = 1:quantity_points
			lattice_distribution(horizontal_points(point), vertical_points(point), 9) = initial_disturbance_density/9;
		end
	end

	lattice{1} = lattice_distribution;		
	
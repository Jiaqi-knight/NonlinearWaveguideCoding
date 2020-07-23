function lattice = set_initial_disturbances(lattice, initial_disturbance_density)
	
	lattice_distribution = lattice{1};
	size_lattice = size(lattice_distribution(:, :, 1));
	number_lines_lattice = size_lattice(1);
	number_columns_lattice = size_lattice(2);
	lattice_distribution(number_lines_lattice/2, number_columns_lattice/2,:) = initial_disturbance_density;
	lattice{1} = lattice_distribution;
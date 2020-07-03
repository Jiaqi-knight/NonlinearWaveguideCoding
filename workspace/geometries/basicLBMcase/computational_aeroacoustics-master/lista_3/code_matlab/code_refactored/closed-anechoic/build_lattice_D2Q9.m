%% functionname: function description
function lattice = build_lattice_D2Q9(number_lines_lattice, number_columns_lattice, lattice_average_density)

% 1
number_directions_D2Q9 = 9;
lattice_distribution = zeros(number_lines_lattice, number_columns_lattice, number_directions_D2Q9);                                 

% 2 - Filling the initial distribution function (at t=0) with initial values
lattice_distribution(:,:,:) = lattice_average_density/9;   
lattice_velocity_x = zeros(number_lines_lattice, number_columns_lattice);
lattice_velocity_y = zeros(number_lines_lattice, number_columns_lattice);

% 3
lattice{1} = lattice_distribution;
lattice{2} = lattice_velocity_x;
lattice{3} = lattice_velocity_y;
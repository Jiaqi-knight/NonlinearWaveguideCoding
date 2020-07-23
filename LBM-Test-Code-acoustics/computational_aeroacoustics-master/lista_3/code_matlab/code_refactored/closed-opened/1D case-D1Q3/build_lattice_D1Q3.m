%% functionname: function description
function lattice = build_lattice_D1Q3( number_columns_lattice, lattice_average_density)

% 1
number_directions_D1Q3= 3;
lattice_distribution = zeros(number_columns_lattice, number_directions_D1Q3);                                 

% 2 - Filling the initial distribution function (at t=0) with initial values
lattice_distribution(:,:,:) = lattice_average_density/3;   
lattice_velocity_x = zeros(1, number_columns_lattice);
lattice_velocity_y = zeros(1, number_columns_lattice);

% 3
lattice{1} = lattice_distribution;
lattice{2} = lattice_velocity_x;
% lattice{3} = lattice_velocity_y;
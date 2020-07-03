% 2D Lattice Boltzmann (BGK) model of a fluid.
%  c6  c2   c5  D2Q9 model. At each timestep, particle densities propagate
%    \  |  /    outwards in the directions indicated in the figure. An
%  c3 -c9 - c1  equivalent 'equilibrium' density is found, and the densities
%    /  |  \    relax towards that state, in a proportion governed by omega.
%  c7  c4   c8      Iain Haslam, March 2006.

clear all, clc
close all

%% 1 - Set lattice sizes
number_lines_lattice = 300; % cells in the y direction
number_columns_lattice = 300; % cells in the x direction

%% 2 - Set physical parameters (macro)
physical_sound_velocity = 340; % [m/s]
physical_density = 1.2; % [kg/m^3]
physical_dimension_max_x = .5; % [m]
physical_dimension_max_y = .5; % [m]
% voxel is a term to express a volume decribed in a pixel: volume + pixel = voxel
dimension_x_voxel = physical_dimension_max_x/number_columns_lattice; % defining dimension x in voxel
lattice_time_step = (1/sqrt(3))*dimension_x_voxel/physical_sound_velocity;

%% 3 - Set lattice parameters (meso - lattice unities)
frequency_relaxation = 1.9; % to 1.5e-5 physcosity 1.9998; 860e-5 = 1.9
time_relaxation = 1/frequency_relaxation;
lattice_average_density = 1;
lattice_sound_speed = 1/sqrt(3);
lattice_sound_speed_pow_2 = lattice_sound_speed^2;
lattice_viscosity = lattice_sound_speed_pow_2*(1/frequency_relaxation-0.5);
physical_viscosity = lattice_viscosity*(dimension_x_voxel^2)/lattice_time_step; % [m^2/s]

% 4 - Build lattice struct with D2Q9
lattice = build_lattice_D2Q9(number_lines_lattice, number_columns_lattice, lattice_average_density);

% 5 - Set initial disturbance
%initial_disturbance_density = 0.01;
%lattice = set_initial_disturbances(lattice, initial_disturbance_density);

%% 6 - Begin the iteractive process
wavelength_discretizations = [8 16 25];
frequency_source = physical_sound_velocity/ ... 
(dimension_x_voxel*wavelength_discretizations(2)); % Hz
amplitude_source = 0.001;
for ta = 1 : 150*sqrt(3)
    
    %% 6.1 - Propagation (streaming)
    lattice = stream_lattice(lattice);

    %% 6.2 - Get density
    lattice_distribution = lattice{1};
    density_total = sum(lattice_distribution,3);
    % Get pressure field in ta = NTS along
    lattice_pressure = 0;
    if ta == 259
    	lattice_pressure = lattice_sound_speed^2*(density_total(150:300, 150) - 1);
    	figure;
		plot(lattice_pressure);
		xlabel('Distance in cells number','FontSize',20);
		ylabel('Lattice Pressure','FontSize',20);
		title('Waves with discretization 16 cells per wavelength','FontSize',20);
		%save lattice_pressure_d_8.mat lattice_pressure;
		phase_wave = 2*pi*frequency_source*(ta - 5)*lattice_time_step
		%physical_viscosity = 1.5e-5;
		%0.001 => pascal
		%[p pos]=cylin_wave();
		%(1/sqrt(3))/20,physical_viscosity,1/sqrt(3),0.001/20,1:150,pi/2
		%[analytical_pressure x] = cylin_wave(frequency_source, physical_viscosity, ...
		%physical_sound_velocity, amplitude_source, 1:number_columns_lattice/2, phase_wave);
		frequency_source_analytical = lattice_sound_speed/ ... 
		(wavelength_discretizations(2));
		[analytical_pressure x] = cylin_wave(frequency_source_analytical, ... 
		physical_viscosity, lattice_sound_speed, ...
		amplitude_source/wavelength_discretizations(2), 1:number_columns_lattice/2, pi/2 + pi/4);

		figure;
		plot(lattice_pressure);
		hold on;
		plot(analytical_pressure, 'r');
		xlabel('Distance in cells number','FontSize',20);
		ylabel('Lattice Pressure','FontSize',20);
		title('Waves with discretization 16 cells per wavelength','FontSize',20);
    end

    %% 6.2.1 - Source sound
    attice_distribution = lattice{1};
	size_lattice = size(lattice_distribution(:, :, 1));
	number_lines_lattice = size_lattice(1);
	number_columns_lattice = size_lattice(2);
    source_sound = lattice_distribution(number_lines_lattice/2, number_columns_lattice/2,9) + ...
    amplitude_source*sin(2*pi*frequency_source*(ta - 1)*lattice_time_step);
	lattice_distribution(number_lines_lattice/2, number_columns_lattice/2,9) = source_sound;
 	lattice = set_initial_disturbances(lattice, source_sound);

    %% 6.3 - Collide
    lattice = collide_lattice(lattice, frequency_relaxation);
    
% Ploting the results in real time   
%surf(density_total - 1), view(2), shading flat, axis equal, caxis([-.00001 .00001])
%grid off
%pause(.00001)
ta
end %  End main time Evolution Loop
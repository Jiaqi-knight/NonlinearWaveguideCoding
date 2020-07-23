clear all;
clc;
close all;

%% 1 - Set lattice sizes
number_lines_lattice = 60 + 2; % cells in the y direction
number_columns_lattice = 250 + 2; % cells in the x direction

%% 2 - Set physical parameters (macro)
physical_sound_velocity = 340; % [m/s]
physical_density = 1.2; % [kg/m^3]
physical_dimension_max_x = 0.3; % [m]
physical_dimension_max_y = 0.072; % [m]
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

% 4 - Build lattice struct with D2Q9 (lattice = Y x X)
lattice = build_lattice_D2Q9(number_lines_lattice, number_columns_lattice, lattice_average_density);

% 4.0.1 - Adding conditions anechoic
distance = 30;
growth_delta = 1;
[sigma_mat9 Ft] = build_anechoic_condition(number_lines_lattice, ... 
number_columns_lattice, distance, growth_delta);

% 4.1 - Setting conditions of wall
wall_points{1} = [2 2 251 251 2]; % set points in horizontal
wall_points{2} = [2 61 61 2 2]; % set points in vertical
conditions_wall = crossing3( number_lines_lattice, ...
number_columns_lattice, wall_points);

% 5 - Set initial disturbance
initial_disturbance_density = 0.001;
points_lattice{2} = [3:60]; % set point along y
points_lattice{1} = linspace(3,3,length(points_lattice{2})); % set all the points in x = 3
lattice = set_initial_disturbances(lattice, initial_disturbance_density, points_lattice);

%% 6 - Begin the iteractive process
time_final = round(4*2*250*sqrt(3))
pressure_final_tube(1:time_final) = 0;
particule_velocity_final_tube(1:time_final) = 0;
for ta = 1 : time_final

    %% 6.1 - Propagation (streaming)
    lattice = stream_lattice(lattice);

    %% 6.1.2 - Setting conditions of walls
    lattice = set_conditions_wall(lattice, conditions_wall);

    %% 6.2 - Recalculating density and velocities
    density = sum(lattice{1},3);
    pressures_input(ta) = mean(density(26:31, 60) - 1)*lattice_sound_speed^2;
    pressures_output(ta) = mean(density(26:31, number_columns_lattice - 60) - 1)*lattice_sound_speed^2;
    lattice = calculate_velocities(lattice, density);

    %% 6.3 - Collide
    lattice = collide_lattice(lattice, frequency_relaxation, ... 
    sigma_mat9, Ft);

    % % Ploting the results in real time
    pressure_final_tube(ta) = mean(density(3:60,250) - 1)*lattice_sound_speed^2;
    horizontal_velocity = lattice{2};
    particule_velocity_final_tube(ta) = mean(horizontal_velocity(3:60,250));  
    %grid off
    imshow(mat2gray(density - 1));
    %imagesc(density-1)
    pause(.0000001) 
    (ta/time_final)*100
end %  End main time Evolution Loop

frequency_pressure_final_tube = fft(pressure_final_tube);
frequency_particule_velocity_final_tube = fft(particule_velocity_final_tube);
impedance = frequency_pressure_final_tube./frequency_particule_velocity_final_tube;
frequencies = linspace(0,1/lattice_time_step,length(frequency_pressure_final_tube))*2*pi*physical_dimension_max_y/physical_sound_velocity;
figure(2);
plot(frequencies, imag(impedance));
ylabel('Impedancia','FontSize',20);
xlabel('Numero de Helmholtz','FontSize',20);
title('Parte Imaginaria da Impedancia','FontSize',20);
axis([0 frequencies(end) -20 20]);
figure(3);
plot(frequencies, real(impedance), 'r');
ylabel('Impedancia','FontSize',20);
xlabel('Numero de Helmholtz','FontSize',20);
title('Parte Real da Impedancia','FontSize',20);
axis([0 frequencies(end) -20 20]);
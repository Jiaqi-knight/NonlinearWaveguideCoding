function ammonite_demo

%Spiral types
spiral_type = 'Logarithmic';
spiral_turns = 5;
points_per_turn = 200;
cross_section_ratio = 1;
helicity = 50;

%Add ridges along ammonite (bumps on elliptical cross section)
add_ridges = 1;
ridge_frequency = 10;
bump_amplitude = 0.01;

%Add bumps along spiral
add_bumps = 0;
spiral_bump_amplitude = 0.5;
spiral_bump_frequency = 10;

%Colouring options
add_ridges_to_colour = 1;
add_bumps_to_colour = 1;

%

%Plot ammonite
figure('name','ammonite demo','color',[ 1 1 1],'renderer','opengl' );
surface_handle = make_ammonite( add_ridges, add_ridges_to_colour, add_bumps,...
    add_bumps_to_colour, spiral_type, spiral_turns, points_per_turn, cross_section_ratio,...
    bump_amplitude, ridge_frequency, spiral_bump_amplitude, spiral_bump_frequency, helicity );
shading interp
axis equal
axis vis3d
camlight;
lighting phong

%Spiral plot
fig = spiral_plot( spiral_type );

%End of code
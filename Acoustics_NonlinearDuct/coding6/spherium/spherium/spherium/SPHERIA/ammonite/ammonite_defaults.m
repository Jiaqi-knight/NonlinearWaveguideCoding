%ammonite_defaults
% Defines default parameters for the ammonite spheria. If S.SPHERIA does
% not have an ammonite field then it will obtain them from this function.

function d = ammonite_defaults

%Spiral types
d.spiral_type = 'Logarithmic';
d.spiral_types = {
    'Logarithmic',...
    'Linear',...
    'Archimedian'
    };
d.spiral_turns = 5;
d.points_per_turn = 200;
d.cross_section_ratio = 0.9;
d.helicity = 0;

%Plot spiral function in a separate window
d.plot_spiral = 0;

%Add ridges along ammonite (bumps on elliptical cross section)
d.add_ridges = 1;
d.ridge_frequency = 10;
d.bump_amplitude = 0.01;

%Add bumps along spiral
d.add_bumps = 0;
d.spiral_bump_amplitude = 0.5;
d.spiral_bump_frequency = 10;

%Colouring options
d.add_ridges_to_colour = 1;
d.add_bumps_to_colour = 1;

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];

%End of code
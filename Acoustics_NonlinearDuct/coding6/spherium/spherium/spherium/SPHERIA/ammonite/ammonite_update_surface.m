%ammonite_update_surface
% Draws a shaded ammonite 3D surface!

function ammonite_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get ammonite data and create new surface
s = m.S.SPHERIA.ammonite;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = make_ammonite( s.add_ridges, s.add_ridges_to_colour, s.add_bumps,...
    s.add_bumps_to_colour, s.spiral_type, s.spiral_turns, s.points_per_turn, s.cross_section_ratio,...
    s.bump_amplitude, s.ridge_frequency, s.spiral_bump_amplitude, s.spiral_bump_frequency, s.helicity );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
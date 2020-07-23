%polyhedron_update_surface
% Draws a shaded polyhedron 3D surface!

function polyhedron_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get polyhedron data and create new surface
s = m.S.SPHERIA.polyhedron;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = make_polyhedron( s.N );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
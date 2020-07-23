%spheria_update_surface
% Draws a shaded spheria 3D surface!

function spheria_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get spheria data and create new surface
s = m.S.SPHERIA.spheria;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = sphericalplot( s.sphere_or_surface, s.spherefunction,  s.N );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
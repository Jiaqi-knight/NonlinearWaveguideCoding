%polyspike_update_surface
% Draws a shaded polyspike 3D surface!

function polyspike_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get polyspike data and create new surface
s = m.S.SPHERIA.polyspike;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = make_polyspike( s.N,s.M,s.k,s.P );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
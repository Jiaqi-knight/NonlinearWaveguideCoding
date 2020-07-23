%klein_update_surface
% Draws a shaded klein 3D surface!

function klein_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get klein data and create new surface
s = m.S.SPHERIA.klein;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = make_klein_bottle_plot( s.pipe_N, s.rot_N, s.rtb, s.rb2sp, s.hcb2sp );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
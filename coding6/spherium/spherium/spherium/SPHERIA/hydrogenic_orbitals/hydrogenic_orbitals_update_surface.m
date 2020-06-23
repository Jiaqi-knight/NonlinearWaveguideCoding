%hydrogenic_orbitals_update_surface
% Draws a shaded hydrogenic_orbital 3D surface!

function hydrogenic_orbitals_update_surface

%Obtain handle to spherium GUI main function and axes
mainfig = findobj('tag','FIGUREspherium');
m = guidata( mainfig );

%Get hydrogenic_orbitals data and create new surface
s = m.S.SPHERIA.hydrogenic_orbitals;
figure(m.FIGUREspherium);
axes(m.AXESspherium);
m.H.surface_handle = make_hydrogenic_orbitals( orbital2L(s.quantum_number_L), s.quantum_number_M, s.N );

%Update spherium main GUI with new surface handle
guidata( mainfig, m );

%End of code
%polyhedron_defaults
% Defines default parameters for the polyhedron spheria. If S.SPHERIA does
% not have an polyhedron field then it will obtain them from this function.

function d = polyhedron_defaults

%Number of vertices in a cross section is 2(N+1)
d.N = 5;

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];

%End of code
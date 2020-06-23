%spheria_defaults
% Defines default parameters for the klein spheria. If S.SPHERIA does
% not have an spheria field then it will obtain them from this function.

function d = spheria_defaults

d.spherefunction = 'abs( cos(elev*azi)*sin(azi*elev) )' ;
d.sphere_or_surfaces = {
    'Sphere',...
    'Surface' };
d.sphere_or_surface = 'Surface';
d.N = 300;

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];

%End of code
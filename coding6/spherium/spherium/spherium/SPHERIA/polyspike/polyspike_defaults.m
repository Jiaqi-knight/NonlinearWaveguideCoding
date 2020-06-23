%polyspike_defaults
% Defines default parameters for the polyspike spheria. If S.SPHERIA does
% not have an polyspike field then it will obtain them from this function.

function d = polyspike_defaults

%Number of azimuth spikes
d.M = 2;

%Number of elevation spikes
d.N = 5;

%Spikiness
d.k = 0;

%Number of parabola points between spikes
d.P = 40;

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];


%End of code
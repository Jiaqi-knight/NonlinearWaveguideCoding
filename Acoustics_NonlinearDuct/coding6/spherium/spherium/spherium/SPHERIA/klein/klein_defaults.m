%klein_defaults
% Defines default parameters for the klein spheria. If S.SPHERIA does
% not have an klein field then it will obtain them from this function.

function d = klein_defaults

% Granularity down the pipe (default 100)
d.pipe_N = 100;

% Rotational granularity (default 50)
d.rot_N = 50;

% Radius of the top bend relative to the small pipe (default 3)
d.rtb = 3;

% Radius of the base relative to the small pipe (default 5)
d.rb2sp = 7;

% Height of the conical base relative to the small pipe (default 10)
d.hcb2sp = 8;

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];

%End of code
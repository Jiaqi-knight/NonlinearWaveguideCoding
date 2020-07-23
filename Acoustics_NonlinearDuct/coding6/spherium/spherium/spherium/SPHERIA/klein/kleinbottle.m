%kleinbottle
% Function which constructs a Klein bottle.
% 
% The (variable) parameters are:
% granularity down the pipe (default 100)
% rotational granularity (default 50)
% radius of the top bend relative to the small pipe (default 3)
% radius of the base relative to the small pipe (default 5)
% height of the conical base relative to the small pipe (default 10)
% 
% LAST UPDATED by Andy French (adapted from klein.m by David Smith,
% obtained from the MATLAB file exchange March 2011).

function [x, y, z] = kleinbottle(varargin)

%Get inputs
l = length(varargin);
switch l
    case 0
        resl = 100;
        resa = 50;
        r2 = 3;
        r3 = 5;
        h1 = 10;
    case 1
        resl = varargin{1};
        resa = 50;
        r2 = 3;
        r3 = 5;
        h1 = 10;
    case 2
        resl = varargin{1};
        resa = varargin{2};
        r2 = 3;
        r3 = 5;
        h1 = 10;
    case 3
        resl = varargin{1};
        resa = varargin{2};
        r2 = varargin{3};
        r3 = 5;
        h1 = 10;
    case 4
        resl = varargin{1};
        resa = varargin{2};
        r2 = varargin{3};
        r3 = 5;
        h1 = 10;
    otherwise
        resl = varargin{1};
        resa = varargin{2};
        r2 = varargin{3};
        r3 = varargin{4};
        h1 = varargin{5};
end

%Define basic grids
ll = linspace(0, 12, resl);
thl = linspace(0, 2*pi, resa);
[l, th] = meshgrid(ll, thl);
bound1 = 3;
bound2 = 6;
bound3 = 9;
bound4 = 12;
s = size(l);
x1 = [];
y1 = [];
z1 = [];
x2 = [];
y2 = [];
z2 = [];
x3 = [];
y3 = [];
z3 = [];
x4 = [];
y4 = [];
z4 = [];

%Define klein bottle parameters

r1 = 1;    % radius of small tube
r4 = (r3 - r1)/2; % radius of lower torus

% first segment
i1 = find(l <= bound1); % indices of first path
l1 = l(i1); % l values
th1 = th(i1); % r values
ph = pi * l1 / r2;  % angle for tube curve
x1 = (cos(ph) - 1) * r2 - r1 * cos(th1);
y1 = r1 * sin(th1);
z1 = l1 * h1 / r2;

%second segment
i2 = find((l > bound1) & (l <= bound2)); % indices of second path
l2 = l(i2) - bound1; % l values
th2 = th(i2); % r values
ph = pi * l2 / r2;  % angle for top curve
x2 = -(r2 + r1*cos(th2)).*cos(ph) - r2;
y2 = r1 * sin(th2);
z2 = (r2 + r1*cos(th2)) .* sin(ph) + h1;

%third segment
i3 = find((l > bound2) & (l <= bound3)); % indices of third path
l3 = l(i3) - bound2; % l values
th3 = th(i3); % r values
ph = pi * l3 / bound1;  % angle for body radius curve
R = r1 + r4*(1 - cos(ph));
x3 = R .* cos(th3);
y3 = R .* sin(th3);
z3 = h1 - l3 * h1 / bound1;

%fourth segment
i4 = find((l > bound3) & (l <= bound4)); % indices of last path
l4 = l(i4) - bound3; % l values
th4 = th(i4); % r values
ph = pi * l4 / bound1;  % angle for closure curve
R = r1 + r4*(1 + cos(ph));
x4 = R .* cos(th4);
y4 = R .* sin(th4);
z4 = -r4 .* sin(ph);

%Construct x,y,z surface
x = reshape([x1' x2' x3' x4']', s);
y = reshape([y1' y2' y3' y4']', s);
z = reshape([z1' z2' z3' z4']', s);

%End of code

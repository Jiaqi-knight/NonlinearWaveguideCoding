%colfromcolormap
% Obtain colours corresponding to value x from curent colormap.

function [R,G,B] = colfromcolormap(x)

%Get vaiable range of current colormap
[xmin,xmax] = caxis;

%Get current colormap
map = colormap;
dim = size(map);

%Make sure x is within the caxis range
x(x>xmax) = xmax;
x(x<xmin) = xmin;

%Interpolate R,G,B values at x
R = interp1( linspace(xmin,xmax,dim(1)),map(:,1).',x );
G = interp1( linspace(xmin,xmax,dim(1)),map(:,2).',x );
B = interp1( linspace(xmin,xmax,dim(1)),map(:,3).',x );

%End of code
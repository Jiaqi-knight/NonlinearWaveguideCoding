%spiral_plot

function fig = spiral_plot( spiral_type )

%%%%%%%%%%%%%%%%%%

%Number of spirals
N = 5;

%%%%%%%%%%%%%%%%%%

% Define spiral function. r is range, t is angle measured clockwise from y
% axis (in radians)
if strcmp( spiral_type, 'Archimedian' )
    spiral_func_str = 'r = (t/2*pi).^2';
elseif strcmp( spiral_type, 'Linear'  )
    spiral_func_str = 'r = t/(2*pi)';
else
    %Logarithmic spiral
    spiral_func_str = 'r = exp(t/(2*pi)) - 1';
end

fig = figure('color',[1 1 1],'name','ammonoite spiral');
t = linspace( 0,2*pi*N,500 );
r = spiral(t,spiral_func_str);
x = r.*sin(t);
y = r.*cos(t);
plot(x,y,'b');
xlabel('x')
ylabel('y')
title(['Spiral function ',spiral_func_str])
axis equal

%%

%Function which evaluates spiral radius from angle t /radians.
function r = spiral(t,spiral_func_str)
eval([ spiral_func_str,';'] )

%End of code
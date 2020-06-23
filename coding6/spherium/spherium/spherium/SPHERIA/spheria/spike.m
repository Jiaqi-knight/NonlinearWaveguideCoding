%spike
% Periodic function which produces an array of spikes.

function y = spike(x,a,b)
x = sin(x);
y = 1 + a*exp(-b*abs(x));

%End of code
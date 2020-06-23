%detonate
% Periodic detotation!

function y = detonate(x,a)
y = 1./abs(sin(x));
y(y>a) = NaN;

%End of code
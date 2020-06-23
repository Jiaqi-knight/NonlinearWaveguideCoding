%make_klein_bottle_plot
% Demonstrates how to plot a Klein bottle.
%
% LAST UPDATED by Andy French. November 2011.

function surface_handle = make_klein_bottle_plot( pipe_N, rot_N, rtb, rb2sp, hcb2sp )

%Contruct Klein bottle, colour proportional to radius
[ x,y,z ] = kleinbottle( pipe_N, rot_N, rtb, rb2sp, hcb2sp );
P = sqrt( x.^2 + y.^2 + z.^2 );
norm = max(max(P));
P = P/norm;
x = x/norm;
y = y/norm;
z = z/norm;
surface_handle=surf(x,y,z,P);
axis equal
axis vis3d
axis tight
axis off
shading interp

%End of code.

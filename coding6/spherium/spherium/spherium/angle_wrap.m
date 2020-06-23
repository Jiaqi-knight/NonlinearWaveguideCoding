%angle_wrap
% Function that wraps an angle x_deg in degrees to the range -180..180
%
% LAST UPDATED by Andy French 14th March 2011

function y_deg = angle_wrap(x_deg)

%Determine remainder by division by 360 deg
y_deg = rem( x_deg, 360 );

%Transform to range -180....180
y_deg( y_deg > 180 ) = y_deg( y_deg > 180 ) - 360;
y_deg( y_deg < -180 ) = y_deg( y_deg < -180 ) + 360;

%%

function angle_wrap_demo

x_deg = linspace(-720,360,1000);
y_deg = angle_wrap(x_deg);
figure; plot(x_deg,y_deg);
xlabel('Angle /deg')
ylabel('Wrapped angle /deg')
title('angle wrapping demo')

%End of code
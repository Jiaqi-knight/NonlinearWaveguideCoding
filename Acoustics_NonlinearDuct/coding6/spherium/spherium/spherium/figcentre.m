%figcenter
% Centres a figure with handle fig in the current display.

function figcentre(fig)

%Get screen dimensions in pixels.
set(0,'units','pixels') ;
screen_size=get(0,'Screensize');
display_size=get(0,'MonitorPositions');

%Get original figure units
orig_units = get(fig,'units');

%Get position vector for figure in pixels
set(fig,'units','pixels');
fig_position=get(gcf,'position') ;

%Define coordinates or figure bottom hand corner from screen bottom hand corner
%that will place the figure window in the center of the screen.
centre_x = 0.5*( display_size(3)-fig_position(3) ) - screen_size(1);
centre_y = 0.5*( display_size(4)-fig_position(4) ) - screen_size(2);

%Modify figure position to centre
set(fig,'position',[centre_x,centre_y,fig_position(3),fig_position(4)] );

%Reset figure units
set(fig,'units',orig_units);

%End of code



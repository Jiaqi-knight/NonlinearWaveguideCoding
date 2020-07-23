%gui_snap
% Function which positions GUI2 underneath GUI1, aligned left. Handles to
% GUI figures and g1 and g2 respectively.

function gui_snap( g1, g2)

%Get outer position coordinates of the two GUI figures, in pixels
set(g1, 'units', 'pixels');
set(g2, 'units', 'pixels');
p1 = get(g1,'outerposition');
p2 = get(g2,'outerposition');

%Move GUI2 underneath GUI1 and align left
p2 = [p1(1),p1(2)-p2(4),p2(3),p2(4)];

%Work out centre of figures
cx = p2(1) + 0.5*max( [p1(3),p2(3)] );
cy = p2(2) + 0.5*( p1(4) + p2(4) );

%Compute deviation of figure centre from centre of display
p0=get(0,'MonitorPositions');
dx = 0.5*p0(3) - cx;
dy = 0.5*p0(4) - cy;

%Shift figures so they are at the centre of the display
set(g1,'outerposition',p1+[dx,dy,0,0])
set(g2,'outerposition',p2+[dx,dy,0,0])

%End of code
%polyhedron_update_gui_from_S
% Updates the polyhedron.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function polyhedron_update_gui_from_S

%Get handle to main spherium GUI figure and polyhedron GUI figure
mainfig = findobj('tag','FIGUREspherium');
polyhedronfig = findobj('tag','FIGUREpolyhedron');
m = guidata( mainfig );
a = guidata( polyhedronfig );

%Get polyhedron GUI settings
if isfield( m.S.SPHERIA,'polyhedron' ) == 0
    %Use defaults
    S.SPHERIA.polyhedron = polyhedron_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.polyhedron;

%Update polyhedron.m GUI
set( a.EDITN,'string', num2str( s.N ) );

%Update gui data
guidata( polyhedronfig, a );

%End of code
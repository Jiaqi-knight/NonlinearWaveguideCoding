%polyhedron_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% polyhedron gui is modified.

function polyhedron_update_S_from_gui

%Get handles to main GUI and polyhedron GUI figures
mainfig = findobj('tag','FIGUREspherium');
polyhedronfig = findobj('tag','FIGUREpolyhedron');
m = guidata( mainfig );
a = guidata( polyhedronfig );

%Update polyhedron parameters within data structure S
m.S.SPHERIA.polyhedron.N = str2num( get( a.EDITN,'string' ) );
 
%Update spherium main gui data
guidata( mainfig, m );

%End of code
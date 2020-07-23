%spheria_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% spheria gui is modified.

function spheria_update_S_from_gui

%Get handles to main GUI and spheria GUI figures
mainfig = findobj('tag','FIGUREspherium');
spheriafig = findobj('tag','FIGUREspheria');
m = guidata( mainfig );
a = guidata( spheriafig );

%Update spheria parameters within data structure S
m.S.SPHERIA.spheria.sphere_or_surfaces = get( a.POPUPMENUsphereorsurface, 'string' );
m.S.SPHERIA.spheria.sphere_or_surface =...
    m.S.SPHERIA.spheria.sphere_or_surfaces{ get(a.POPUPMENUsphereorsurface,'value') };
m.S.SPHERIA.spheria.spherefunction = get( a.EDITspherefunction, 'string' );
m.S.SPHERIA.spheria.N = str2num( get( a.EDITN, 'string' ) );
 
%Update spherium main gui data
guidata( mainfig, m );

%End of code
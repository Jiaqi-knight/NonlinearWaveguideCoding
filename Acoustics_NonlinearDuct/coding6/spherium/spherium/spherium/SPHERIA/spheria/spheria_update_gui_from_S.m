%spheria_update_gui_from_S
% Updates the spheria.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function spheria_update_gui_from_S

%Get handle to main spherium GUI figure and spheria GUI figure
mainfig = findobj('tag','FIGUREspherium');
spheriafig = findobj('tag','FIGUREspheria');
m = guidata( mainfig );
a = guidata( spheriafig );

%Get spheria GUI settings
if isfield( m.S.SPHERIA,'spheria' ) == 0
    %Use defaults
    S.SPHERIA.spheria = spheria_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.spheria;

%Update spheria.m GUI
set( a.POPUPMENUsphereorsurface, 'string', s.sphere_or_surfaces );
set( a.POPUPMENUsphereorsurface,'value',...
    strmatch( s.sphere_or_surface, s.sphere_or_surfaces, 'exact' ));
set( a.EDITspherefunction,'string',s.spherefunction );
set( a.EDITN,'string',num2str( s.N ) );

%Update gui data
guidata( spheriafig, a );

%End of code
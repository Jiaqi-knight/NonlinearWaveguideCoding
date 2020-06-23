%polyspike_update_gui_from_S
% Updates the polyspike.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function polyspike_update_gui_from_S

%Get handle to main spherium GUI figure and polyspike GUI figure
mainfig = findobj('tag','FIGUREspherium');
polyspikefig = findobj('tag','FIGUREpolyspike');
m = guidata( mainfig );
a = guidata( polyspikefig );

%Get polyspike GUI settings
if isfield( m.S.SPHERIA,'polyspike' ) == 0
    %Use defaults
    S.SPHERIA.polyspike = polyspike_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.polyspike;

%Update polyspike.m GUI
set( a.EDITaziM,'string', num2str( s.M ) );
set( a.EDITelevN,'string', num2str( s.N ) );
set( a.EDITk,'string', num2str( s.k ) );
set( a.EDITP,'string', num2str( s.P ) );

%Update gui data
guidata( polyspikefig, a );

%End of code
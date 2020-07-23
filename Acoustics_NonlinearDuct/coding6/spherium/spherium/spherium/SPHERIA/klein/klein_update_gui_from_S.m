%klein_update_gui_from_S
% Updates the klein.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function klein_update_gui_from_S

%Get handle to main spherium GUI figure and klein GUI figure
mainfig = findobj('tag','FIGUREspherium');
kleinfig = findobj('tag','FIGUREklein');
m = guidata( mainfig );
a = guidata( kleinfig );

%Get klein GUI settings
if isfield( m.S.SPHERIA,'klein' ) == 0
    %Use defaults
    S.SPHERIA.klein = klein_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.klein;

%Update klein.m GUI
set( a.EDITpipeN,'string', num2str( s.pipe_N ) );
set( a.EDITrotN,'string', num2str( s.rot_N ) );
set( a.EDITrtb,'string', num2str( s.rtb ) );
set( a.EDITrb2sp,'string', num2str( s.rb2sp ) );
set( a.EDIThcb2sp,'string', num2str( s.hcb2sp ) );

%Update gui data
guidata( kleinfig, a );

%End of code
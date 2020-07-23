%klein_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% klein gui is modified.

function klein_update_S_from_gui

%Get handles to main GUI and klein GUI figures
mainfig = findobj('tag','FIGUREspherium');
kleinfig = findobj('tag','FIGUREklein');
m = guidata( mainfig );
a = guidata( kleinfig );

%Update klein parameters within data structure S
m.S.SPHERIA.klein.pipe_N = str2num( get( a.EDITpipeN,'string' ) );
m.S.SPHERIA.klein.rot_N = str2num( get( a.EDITrotN,'string' ) );
m.S.SPHERIA.klein.rtb = str2num( get( a.EDITrtb,'string' ) );
m.S.SPHERIA.klein.rb2sp = str2num( get( a.EDITrb2sp,'string' ) );
m.S.SPHERIA.klein.hcb2sp = str2num( get( a.EDIThcb2sp,'string' ) );
 
%Update spherium main gui data
guidata( mainfig, m );

%End of code
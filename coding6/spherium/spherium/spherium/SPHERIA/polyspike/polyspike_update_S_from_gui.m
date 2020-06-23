%polyspike_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% polyspike gui is modified.

function polyspike_update_S_from_gui

%Get handles to main GUI and polyspike GUI figures
mainfig = findobj('tag','FIGUREspherium');
polyspikefig = findobj('tag','FIGUREpolyspike');
m = guidata( mainfig );
a = guidata( polyspikefig );

%Update polyspike parameters within data structure S
m.S.SPHERIA.polyspike.N = str2num( get( a.EDITelevN,'string' ) );
m.S.SPHERIA.polyspike.M = str2num( get( a.EDITaziM,'string' ) );
m.S.SPHERIA.polyspike.k = str2num( get( a.EDITk,'string' ) );
m.S.SPHERIA.polyspike.P = str2num( get( a.EDITP,'string' ) );
 
%Update spherium main gui data
guidata( mainfig, m );

%End of code
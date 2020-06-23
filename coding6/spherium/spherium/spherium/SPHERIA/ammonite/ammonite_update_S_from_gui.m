%ammonite_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% ammonite gui is modified.

function ammonite_update_S_from_gui

%Get handles to main GUI and ammonite GUI figures
mainfig = findobj('tag','FIGUREspherium');
ammonitefig = findobj('tag','FIGUREammonite');
m = guidata( mainfig );
a = guidata( ammonitefig );

%Update ammonite parameters within data structure S
m.S.SPHERIA.ammonite.spiral_types = get( a.POPUPMENUspiraltype, 'string' );
m.S.SPHERIA.ammonite.spiral_type = m.S.SPHERIA.ammonite.spiral_types{ get(a.POPUPMENUspiraltype,'value') };
m.S.SPHERIA.ammonite.spiral_turns = str2num( get( a.EDITspiralturns, 'string' ) );
m.S.SPHERIA.ammonite.points_per_turn = str2num( get( a.EDITpointsperturn, 'string' ) );
m.S.SPHERIA.ammonite.cross_section_ratio = str2num( get( a.EDITcrosssectionratio, 'string' ) );
m.S.SPHERIA.ammonite.helicity = str2num( get( a.EDIThelicity, 'string' ) );
m.S.SPHERIA.ammonite.plot_spiral = get( a.CHECKplotspiral,'value' );
m.S.SPHERIA.ammonite.add_ridges = get( a.CHECKaddridges,'value' );
m.S.SPHERIA.ammonite.ridge_frequency = str2num( get( a.EDITridgefrequency, 'string' ) );
m.S.SPHERIA.ammonite.bump_amplitude = str2num( get(a.EDITbumpamplitude,'string' ) );
m.S.SPHERIA.ammonite.add_bumps = get( a.CHECKaddbumps,'value' );
m.S.SPHERIA.ammonite.spiral_bump_frequency = str2num( get( a.EDITspiralbumpfrequency, 'string' ) );
m.S.SPHERIA.ammonite.spiral_bump_amplitude = str2num( get( a.EDITspiralbumpamplitude, 'string' ) );
m.S.SPHERIA.ammonite.add_ridges_to_colour = get( a.CHECKaddridgestocolour,'value' );
m.S.SPHERIA.ammonite.add_bumps_to_colour = get( a.CHECKaddbumpstocolour,'value' );

%Update spherium main gui data
guidata( mainfig, m );

%End of code
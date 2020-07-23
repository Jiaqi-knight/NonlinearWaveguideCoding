%ammonite_update_gui_from_S
% Updates the ammonite.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function ammonite_update_gui_from_S

%Get handle to main spherium GUI figure and ammonite GUI figure
mainfig = findobj('tag','FIGUREspherium');
ammonitefig = findobj('tag','FIGUREammonite');
m = guidata( mainfig );
a = guidata( ammonitefig );

%Get ammonite GUI settings
if isfield( m.S.SPHERIA,'ammonite' ) == 0
    %Use defaults
    S.SPHERIA.ammonite = ammonite_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.ammonite;

%Spiral types
set( a.POPUPMENUspiraltype, 'string', s.spiral_types );
set( a.POPUPMENUspiraltype,'value',...
    strmatch( s.spiral_type, s.spiral_types, 'exact' ));
set( a.EDITspiralturns,'string',num2str(s.spiral_turns) );
set( a.EDITpointsperturn,'string',num2str(s.points_per_turn) );
set( a.EDITcrosssectionratio,'string',num2str(s.cross_section_ratio) );
set( a.EDIThelicity,'string',num2str(s.helicity) );

%Plot spiral function in a separate window
set( a.CHECKplotspiral,'value',s.plot_spiral );

%Add ridges along ammonite (bumps on elliptical cross section)
set( a.CHECKaddridges,'value',s.add_ridges ) ;
set( a.EDITridgefrequency,'string',num2str(s.ridge_frequency) ) ;
set( a.EDITbumpamplitude,'string',num2str(s.bump_amplitude) ) ;

%Add bumps along spiral
set( a.CHECKaddbumps,'value',s.add_bumps );
set( a.EDITspiralbumpamplitude,'string',...
    num2str(s.spiral_bump_amplitude) ) ;
set( a.EDITspiralbumpfrequency,'string',num2str(s.spiral_bump_frequency) ) ;

%Colouring options
set( a.CHECKaddridgestocolour,'value',s.add_ridges_to_colour ) ;
set( a.CHECKaddbumpstocolour,'value',s.add_bumps_to_colour ) ;

%Update gui data
guidata( ammonitefig, a );

%End of code
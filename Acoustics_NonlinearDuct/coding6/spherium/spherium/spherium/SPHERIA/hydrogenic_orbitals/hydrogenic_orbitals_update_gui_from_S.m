%hydrogenic_orbitals_update_gui_from_S
% Updates the hydrogenic_orbitals.m gui based upon parameters contained within data
% structure S held within the gui data of the spherium main gui.

function hydrogenic_orbitals_update_gui_from_S

%Get handle to main spherium GUI figure and hydrogenic_orbitals GUI figure
mainfig = findobj('tag','FIGUREspherium');
hydrogenic_orbitalsfig = findobj('tag','FIGUREhydrogenic_orbitals');
m = guidata( mainfig );
a = guidata( hydrogenic_orbitalsfig );

%Get hydrogenic_orbitals GUI settings
if isfield( m.S.SPHERIA,'hydrogenic_orbitals' ) == 0
    %Use defaults
    S.SPHERIA.hydrogenic_orbitals = hydrogenic_orbitals_defaults;
end

%Copy parameters to structure s for brevity
s = m.S.SPHERIA.hydrogenic_orbitals;

%Update hydrogenic_orbitals.m GUI
set( a.POPUPMENUquantumnumN, 'string', s.quantum_number_Ns );
set( a.POPUPMENUquantumnumN,'value',...
    strmatch( num2str(s.quantum_number_N), s.quantum_number_Ns, 'exact' ));
set( a.POPUPMENUquantumnumL, 'string', s.quantum_number_Ls );
set( a.POPUPMENUquantumnumL,'value',...
    strmatch( s.quantum_number_L, s.quantum_number_Ls, 'exact' ));
set( a.POPUPMENUquantumnumM, 'string', s.quantum_number_Ms );
set( a.POPUPMENUquantumnumM,'value',...
    strmatch( num2str(s.quantum_number_M), s.quantum_number_Ms, 'exact' ));
set( a.EDITZ,'string',num2str( s.Z ) );
set( a.EDITA,'string',num2str( s.A ) );
set( a.EDITenergy,'string',num2str( s.E_eV ) );
set( a.EDITN,'string',num2str( s.N ) );

%Update gui data
guidata( hydrogenic_orbitalsfig, a );

%End of code
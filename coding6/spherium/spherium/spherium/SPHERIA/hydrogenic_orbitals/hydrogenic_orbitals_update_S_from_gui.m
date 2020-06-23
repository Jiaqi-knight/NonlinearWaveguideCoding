%hydrogenic_orbitals_update_S_from_gui
% Updades data structure S within gui data of main spherium gui when the
% hydrogenic_orbitals gui is modified.

function hydrogenic_orbitals_update_S_from_gui

%Get handles to main GUI and hydrogenic_orbitals GUI figures
mainfig = findobj('tag','FIGUREspherium');
hydrogenic_orbitalsfig = findobj('tag','FIGUREhydrogenic_orbitals');
m = guidata( mainfig );
a = guidata( hydrogenic_orbitalsfig );

%Update hydrogenic_orbitals parameters within data structure S
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ns = get( a.POPUPMENUquantumnumN, 'string' );
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_N =...
    str2num( m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ns{ get(a.POPUPMENUquantumnumN,'value') });
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ls = get( a.POPUPMENUquantumnumL, 'string' );
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_L =...
    m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ls{ get(a.POPUPMENUquantumnumL,'value') };
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ms = get( a.POPUPMENUquantumnumM, 'string' );
m.S.SPHERIA.hydrogenic_orbitals.quantum_number_M =...
    str2num( m.S.SPHERIA.hydrogenic_orbitals.quantum_number_Ms{ get(a.POPUPMENUquantumnumM,'value') });
m.S.SPHERIA.hydrogenic_orbitals.Z = str2num( get( a.EDITZ, 'string' ) );
m.S.SPHERIA.hydrogenic_orbitals.A = str2num( get( a.EDITA, 'string' ) );
m.S.SPHERIA.hydrogenic_orbitals.E_eV = str2num( get( a.EDITenergy, 'string' ) );
m.S.SPHERIA.hydrogenic_orbitals.N = str2num( get( a.EDITN, 'string' ) );
 
%Update spherium main gui data
guidata( mainfig, m );

%End of code
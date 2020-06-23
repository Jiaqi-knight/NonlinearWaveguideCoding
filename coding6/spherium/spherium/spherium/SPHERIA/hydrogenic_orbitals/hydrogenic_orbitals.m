%hydrogenic_orbitals
% GUI front end for hydrogenic_orbitals.m

function varargout = hydrogenic_orbitals(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @hydrogenic_orbitals_OpeningFcn, ...
    'gui_OutputFcn',  @hydrogenic_orbitals_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

%%

%% hydrogenic_orbitals opening function. This executes just before hydrogenic_orbitals is made visible.
function hydrogenic_orbitals_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
hydrogenic_orbitals_update_gui_from_S

%%

%% hydrogenic_orbitals output function. Outputs from this function are returned to the command line.
function varargout = hydrogenic_orbitals_OutputFcn(hObject, eventdata, handles)

%%

%% Number of surface points edit box
function EDITN_Callback(hObject, eventdata, handles)
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%% Atomic number edit box
function EDITZ_Callback(hObject, eventdata, handles)

%Update energy calculation
handles = guidata(hObject);
Ns = get( handles.POPUPMENUquantumnumN,'string' );
N = str2num( Ns{get( handles.POPUPMENUquantumnumN,'value' )} );
Ls = get( handles.POPUPMENUquantumnumL,'string' );
L = Ls{get( handles.POPUPMENUquantumnumL,'value' )};
Z = str2num( get( handles.EDITZ,'string' ) );
A = str2num( get( handles.EDITA,'string' ) );
[w,E_eV,r_mean,r2_mean,a0,a] = Hradial(0,N,orbital2L(L),Z,A);

%Update GUI
set( handles.EDITenergy,'string',num2str(E_eV) );
guidata(hObject,handles);
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%% Atomic mass number edit box
function EDITA_Callback(hObject, eventdata, handles)

%Update energy calculation
handles = guidata(hObject);
Ns = get( handles.POPUPMENUquantumnumN,'string' );
N = str2num( Ns{get( handles.POPUPMENUquantumnumN,'value' )} );
Ls = get( handles.POPUPMENUquantumnumL,'string' );
L = Ls{get( handles.POPUPMENUquantumnumL,'value' )};
Z = str2num( get( handles.EDITZ,'string' ) );
A = str2num( get( handles.EDITA,'string' ) );
[w,E_eV,r_mean,r2_mean,a0,a] = Hradial(0,N,orbital2L(L),Z,A);

%Update GUI
set( handles.EDITenergy,'string',num2str(E_eV) );
guidata(hObject,handles);
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%% Quantum number N pop-up menu
function POPUPMENUquantumnumN_Callback(hObject, eventdata, handles)

%Change L and M quantum number values appropriate to N selection
handles = guidata(hObject);
Ns = get( handles.POPUPMENUquantumnumN,'string' );
N = str2num( Ns{get( handles.POPUPMENUquantumnumN,'value' )} );
Ls = LgivenN(N);
L = Ls{1};
Ms = MgivenL(L);

%Update energy calculation
Z = str2num( get( handles.EDITZ,'string' ) );
A = str2num( get( handles.EDITA,'string' ) );
[w,E_eV,r_mean,r2_mean,a0,a] = Hradial(0,N,orbital2L(L),Z,A);

%Update GUI
set( handles.POPUPMENUquantumnumL,'value',1,'string',Ls );
set( handles.POPUPMENUquantumnumM,'value',1,'string',Ms );
set( handles.EDITenergy,'string',num2str(E_eV) );
guidata(hObject,handles);
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%% Orbital quantum number L pop-up menu
function POPUPMENUquantumnumL_Callback(hObject, eventdata, handles)

%Change M quantum number values appropriate to L selection
handles = guidata(hObject);
Ns = get( handles.POPUPMENUquantumnumN,'string' );
N = str2num( Ns{get( handles.POPUPMENUquantumnumN,'value' )} );
Ls = get( handles.POPUPMENUquantumnumL,'string' );
L = Ls{get( handles.POPUPMENUquantumnumL,'value' )};
Ms = MgivenL(L);

%Update energy calculation
Z = str2num( get( handles.EDITZ,'string' ) );
A = str2num( get( handles.EDITA,'string' ) );
[w,E_eV,r_mean,r2_mean,a0,a] = Hradial(0,N,orbital2L(L),Z,A);

%Update GUI
set( handles.POPUPMENUquantumnumM,'value',1,'string',Ms );
set( handles.EDITenergy,'string',num2str(E_eV) );
guidata(hObject,handles);
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%% Magnetic quantum number M pop-up menu
function POPUPMENUquantumnumM_Callback(hObject, eventdata, handles)
hydrogenic_orbitals_update_S_from_gui;
hydrogenic_orbitals_update_surface
spherium_update;

%%

%End of code







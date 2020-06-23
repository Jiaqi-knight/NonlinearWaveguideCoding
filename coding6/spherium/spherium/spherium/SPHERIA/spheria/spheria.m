%spheria
% GUI front end for spheria.m

function varargout = spheria(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spheria_OpeningFcn, ...
    'gui_OutputFcn',  @spheria_OutputFcn, ...
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

%% spheria opening function. This executes just before spheria is made visible.
function spheria_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
spheria_update_gui_from_S

%%

%% spheria output function. Outputs from this function are returned to the command line.
function varargout = spheria_OutputFcn(hObject, eventdata, handles)

%%

%% Number of surface points edit box
function EDITN_Callback(hObject, eventdata, handles)
spheria_update_S_from_gui;
spheria_update_surface
spherium_update;

%%

%% Spherefunction edit box
function EDITspherefunction_Callback(hObject, eventdata, handles)
spheria_update_S_from_gui;
spheria_update_surface
spherium_update;

%%

%% Select surface or sphere surface option
function POPUPMENUsphereorsurface_Callback(hObject, eventdata, handles)
spheria_update_S_from_gui;
spheria_update_surface
spherium_update;

%%

%End of code







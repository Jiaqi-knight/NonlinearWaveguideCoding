%ammonite
% GUI front end for ammonite.m

function varargout = ammonite(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ammonite_OpeningFcn, ...
    'gui_OutputFcn',  @ammonite_OutputFcn, ...
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

%% ammonite opening function. This executes just before ammonite is made visible.
function ammonite_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
ammonite_update_gui_from_S

%%

%% ammonite output function. Outputs from this function are returned to the command line.
function varargout = ammonite_OutputFcn(hObject, eventdata, handles)

%%

%% Ammonite plot spiral checkbox
function CHECKplotspiral_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite add ridges checkbox
function CHECKaddridges_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite add bumps checkbox
function CHECKaddbumps_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite add ridges to colour checkbox
function CHECKaddridgestocolour_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite add bumps to colour checkbox
function CHECKaddbumpstocolour_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite spiral type popupmenu
function POPUPMENUspiraltype_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite # spiral turns edit box
function EDITspiralturns_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite # surface points per turn edit box
function EDITpointsperturn_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite ridge frequency edit box
function EDITridgefrequency_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite cross section ratio edit box
function EDITcrosssectionratio_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite bump amplitude edit box
function EDITbumpamplitude_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite spiral bump amplitude edit box
function EDITspiralbumpamplitude_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite helicity edit box
function EDIThelicity_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%% Ammonite spiral bump frequency edit box
function EDITspiralbumpfrequency_Callback(hObject, eventdata, handles)
ammonite_update_S_from_gui;
ammonite_update_surface
spherium_update;

%%

%End of code







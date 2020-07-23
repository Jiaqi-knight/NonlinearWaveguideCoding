%klein
% GUI front end for klein.m

function varargout = klein(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @klein_OpeningFcn, ...
    'gui_OutputFcn',  @klein_OutputFcn, ...
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

%% klein opening function. This executes just before klein is made visible.
function klein_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
klein_update_gui_from_S

%%

%% klein output function. Outputs from this function are returned to the command line.
function varargout = klein_OutputFcn(hObject, eventdata, handles)

%%

%% Pipe granularity edit box
function EDITpipeN_Callback(hObject, eventdata, handles)
klein_update_S_from_gui;
klein_update_surface
spherium_update;

%%

%% Rotational granulaity edit box
function EDITrotN_Callback(hObject, eventdata, handles)
klein_update_S_from_gui;
klein_update_surface
spherium_update;

%%

%% Radius of top bend to small pipe ratio edit box
function EDITrtb_Callback(hObject, eventdata, handles)
klein_update_S_from_gui;
klein_update_surface
spherium_update;

%%

%% Radius of base to small pipe ratio edit box
function EDITrb2sp_Callback(hObject, eventdata, handles)
klein_update_S_from_gui;
klein_update_surface
spherium_update;

%%

%% Cone height to small pipe radius ratio edit box
function EDIThcb2sp_Callback(hObject, eventdata, handles)
klein_update_S_from_gui;
klein_update_surface
spherium_update;

%%

%End of code







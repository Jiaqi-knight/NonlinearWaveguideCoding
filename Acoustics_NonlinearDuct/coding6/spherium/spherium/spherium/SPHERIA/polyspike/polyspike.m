%polyspike
% GUI front end for polyspike.m

function varargout = polyspike(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @polyspike_OpeningFcn, ...
    'gui_OutputFcn',  @polyspike_OutputFcn, ...
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

%% polyspike opening function. This executes just before polyspike is made visible.
function polyspike_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
polyspike_update_gui_from_S

%%

%% polyspike output function. Outputs from this function are returned to the command line.
function varargout = polyspike_OutputFcn(hObject, eventdata, handles)

%%

%% # azimuth spikes
function EDITaziM_Callback(hObject, eventdata, handles)
polyspike_update_S_from_gui;
polyspike_update_surface
spherium_update;

%%

%% # elevation spikes
function EDITelevN_Callback(hObject, eventdata, handles)
polyspike_update_S_from_gui;
polyspike_update_surface
spherium_update;

%%

%% Spikiness
function EDITk_Callback(hObject, eventdata, handles)
polyspike_update_S_from_gui;
polyspike_update_surface
spherium_update;

%%

%% # points per spike parabola
function EDITP_Callback(hObject, eventdata, handles)
polyspike_update_S_from_gui;
polyspike_update_surface
spherium_update;

%%

%End of code







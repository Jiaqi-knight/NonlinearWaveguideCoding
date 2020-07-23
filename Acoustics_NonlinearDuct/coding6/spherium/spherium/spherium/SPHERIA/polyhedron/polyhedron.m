%polyhedron
% GUI front end for polyhedron.m

function varargout = polyhedron(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @polyhedron_OpeningFcn, ...
    'gui_OutputFcn',  @polyhedron_OutputFcn, ...
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

%% polyhedron opening function. This executes just before polyhedron is made visible.
function polyhedron_OpeningFcn(hObject, eventdata, handles, varargin)

%Set GUI parameters from data structure S held within the GUI data of the
%main spherium GUI.
polyhedron_update_gui_from_S

%%

%% polyhedron output function. Outputs from this function are returned to the command line.
function varargout = polyhedron_OutputFcn(hObject, eventdata, handles)

%%

%% Number of vertices in a cross section is 2(N+1)
function EDITN_Callback(hObject, eventdata, handles)
polyhedron_update_S_from_gui;
polyhedron_update_surface
spherium_update;

%End of code







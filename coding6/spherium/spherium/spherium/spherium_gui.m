%spherium_gui
% GUI front end for spherium.m

function varargout = spherium_gui(varargin)

%% Initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spherium_gui_OpeningFcn, ...
    'gui_OutputFcn',  @spherium_gui_OutputFcn, ...
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

%% spherium opening function. This executes just before spherium is made visible.
function spherium_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for spherium
handles.output = hObject;

%Ensure plot axis starts blank
axes( handles.AXESspherium );
axis off;

%Center figure
figcentre( handles.FIGUREspherium );

%Turn off warnings (e.g. 'divide by zero')
warning off

%Load spherium logo
axes( handles.AXESlogo );
image( imread('spherium\spherium_logo.png' ));
axis tight
axis equal
axis off

% Load previous GUI settings, or default
[handles.S,handles.H] = spherium_defaults;
if exist('spherium\last_spherium_settings.mat','file')
    load('spherium\last_spherium_settings.mat');
    handles.S = S;
    clear S
end
handles.H.colorbar_handle = [];
handles.H.surface_handle = [];

%Load other data in spherium_store and store this in the handles array
load('spherium\spherium_store');
handles.last_lookin_dir = last_lookin_dir;
if exist(last_lookin_dir,'dir')==0
    handles.last_lookin_dir = pwd;
end

%Initialise mouse button status
handles.button = 'not clicked';

%Set GUI flag that it is open and ready to accept windowbuttonmove
%callbacks
handles.GUIisopen = 1;

%Update GUI
guidata( hObject, handles );
spherium_update_gui_from_S;

%Update spherium
spherium_update

%Load current spheria parameter GUI
eval( handles.S.spheria );
h = findobj('tag',['FIGURE',handles.S.spheria]);

%Align guis
gui_snap( handles.FIGUREspherium,h );

%%

%% spherium output function. Outputs from this function are returned to the command line.
function varargout = spherium_gui_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%

%% Save .PNG file of current spherium
function PUSHsavepng_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(hObject);

%Start clock
t0=cputime;

%Update status
handles.S.status = spherium_status('Saving PNG');
set( handles.TEXTstatus,'string',handles.S.status );

%Obtain colormap and colour limits of spherium axes
axes(handles.AXESspherium);
cmap = colormap;
cax = caxis;

%Create new figure containing a copy of the spherium axes. Save a .png
%file bitmap from this figure into and then close the figure.
new_fig = figure('color',[1 1 1],'name','spherium png output',...
    'units','normalized');
new_ax = copyobj(handles.AXESspherium,new_fig);
if ~isempty(handles.H.colorbar_handle)
    set( new_fig, 'position', [0 0 1 1] );
    new_cbar = copyobj(handles.H.colorbar_handle,new_fig);
end
colormap(cmap);
caxis(cax)

%Crop figure to image and centre
if handles.S.add_axis == 0
    cropfig( new_fig, new_ax, 2 );
else
    cropfig( new_fig, new_ax, 1 );
end
figcentre(new_fig);

%Determine image height in cm
set(new_fig,'units','centimeters');
fp = get(new_fig,'position');
fig_height_cm = fp(4);

%Determine output image width and height in cm
[width_cm,height_cm] = papersize( handles.S.output_size );

%Determine DPI such that output image height at handles.S.output_DPI is
%height_cm
DPI = handles.S.output_DPI * height_cm / fig_height_cm ;

%Print date and time labelled image at required DPI
filename = ['spherium image ',strrep(datestr(now),':','-')];
print(new_fig,filename,'-dpng', ['-r',num2str(handles.S.output_DPI)] );

%Crop image to non-white regions
croppic( [filename,'.png'],[255,255,255] );

%Add blanks to image to fit onto desired paper size
fitpic2paper( [filename,'.png'],width_cm/height_cm, [255,255,255] );
close(new_fig);

%Restore GUI with axis, colorbars etc. Also force light toggle to be off.
handles.S.light_toggle = 0;
guidata(hObject,handles);
spherium_update;

%Update status
handles.S.status = ['PNG image saved in ',num2str(cputime-t0),' s.  ',...
    spherium_status('default') ];
set( handles.TEXTstatus,'string',handles.S.status );

%%

%% Load default settings pushbutton
function PUSHdefault_Callback(hObject, eventdata, handles)

%Get spherium default settings
S = spherium_defaults;

%Store spherium settings in handles array, to enable GUI callback
%functions to see this data
handles = guidata(hObject);
handles.S = S;

%Reset camera position and target
set( handles.AXESspherium,'CameraPosition',handles.S.camera_position,...
    'CameraTarget',handles.S.camera_target,'CameraViewAngle',handles.S.camera_view_angle,...
    'CameraUpVector',handles.S.camera_up_vector);

%Get axes properties and store these
[handles.S.axes_properties_fields,handles.S.axes_properties] =...
    getax(handles.AXESspherium);

%Copy default holes, texture and transparency maps
copyfile( 'spherium\default_texture.jpg','spherium\texture.jpg','f');
copyfile( 'spherium\default_holes.jpg','spherium\holes.jpg','f');
copyfile( 'spherium\default_transparency.jpg','spherium\transparency.jpg','f');

%Update GUI
guidata(hObject,handles);
spherium_update_gui_from_S;
spherium_update;

%%

%% Azimuth /deg edit box
function EDITazi_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get azimuth /deg from editbox and convert to range -180...180
handles.S.azi_deg = str2num( get(handles.EDITazi,'string') );
handles.S.azi_deg = angle_wrap( handles.S.azi_deg ) ;

%Rotate spheria
p = get( handles.AXESspherium,'CameraPosition' );
[x,y,z] = sph2cart( handles.S.azi_deg*pi/180,...
    handles.S.elev_deg*pi/180,sqrt(p(1)^2+p(2)^2+p(3)^2) );
set( handles.AXESspherium, 'CameraPosition',[x,y,z] );
guidata(gcf,handles);
spherium_fast_update;

%%

%% Azimuth /deg slider
function SLIDERazi_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value from slider and convert to range -180...180
handles.S.azi_deg = -180 + get(handles.SLIDERazi,'value')*360;
p = get( handles.AXESspherium,'CameraPosition' );
[x,y,z] = sph2cart( handles.S.azi_deg*pi/180,...
    handles.S.elev_deg*pi/180,sqrt(p(1)^2+p(2)^2+p(3)^2) );
set( handles.AXESspherium, 'CameraPosition',[x,y,z] );
guidata(gcf,handles);
spherium_fast_update;

%%

%% Elevation /deg edit box
function EDITelev_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get elev in /deg from editbox and convert to range -90...90
handles.S.elev_deg = str2num( get(handles.EDITelev,'string') );
handles.S.elev_deg = angle_wrap( handles.S.elev_deg ) ;

%Prevent any angles greater than 90 or less than -90
if handles.S.elev_deg > 90
    handles.S.elev_deg = 90;
elseif handles.S.elev_deg < -90
    handles.S.elev_deg = -90;
end

%Update GUI
p = get( handles.AXESspherium,'CameraPosition' );
[x,y,z] = sph2cart( handles.S.azi_deg*pi/180,...
    handles.S.elev_deg*pi/180,sqrt(p(1)^2+p(2)^2+p(3)^2) );
set( handles.AXESspherium, 'CameraPosition',[x,y,z] );
guidata(gcf,handles);
spherium_fast_update;

%%

%% Azimuth /deg slider
function SLIDERelev_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value from slider and convert to range -90...90
handles.S.elev_deg = -90 + get(handles.SLIDERelev,'value')*180;

%Update GUI
p = get( handles.AXESspherium,'CameraPosition' );
[x,y,z] = sph2cart( handles.S.azi_deg*pi/180,...
    handles.S.elev_deg*pi/180,sqrt(p(1)^2+p(2)^2+p(3)^2) );
set( handles.AXESspherium, 'CameraPosition',[x,y,z] );
guidata(gcf,handles);
spherium_fast_update;

%%

%% Output image size popup menu
function POPUPMENUoutputsize_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get string from menu
v = get( handles.POPUPMENUoutputsize, 'value' );
handles.S.output_size = handles.S.output_sizes{v};

%Update GUI
guidata(hObject,handles);

%%

%% Plot type popup menu
function POPUPMENUspheria_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Close previous GUI
delete( findobj('tag',['FIGURE',handles.S.spheria] ));

%Get string from menu
v = get( handles.POPUPMENUspheria, 'value' );
handles.S.spheria = handles.S.spherias{v};

%Load new GUI
eval( handles.S.spheria );
h = findobj('tag',['FIGURE',handles.S.spheria]);

%Align guis
gui_snap( handles.FIGUREspherium,h );

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%% Output image Dots Per Inch edit box
function EDITDPI_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get DPI
handles.S.output_DPI = str2num( get(handles.EDITDPI,'string') );
handles.S.output_DPI = handles.S.output_DPI(1);

%Update GUI
guidata(hObject,handles);

%%

%% Renderer popupmenu
function POPUPMENUrenderer_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get string from menu
v = get( handles.POPUPMENUrenderer, 'value' );
handles.S.renderer = handles.H.renderers{v};

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%% Colour map popupmenu
function POPUPMENUcolormap_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get string from menu
v = get( handles.POPUPMENUcolormap, 'value' );
handles.S.colour_map = handles.H.colour_maps{v};

%Load bespoke colormap if this option is set
if strcmp(handles.S.colour_map,'bespoke')==1
    [filename, pathname, filterindex] = uigetfile('*.mat', 'Select bespoke colormap file');
    if filename == 0
        load( 'bespoke_colormap_default.mat');
    else
        load([pathname,'\',filename]);
    end
    
    %Set colormap
    colormap(map);
end

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

function POPUPMENUcolorfunction_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get string from menu
handles.S.colorfunction = ...
    handles.S.colorfunctions{ get(handles.POPUPMENUcolorfunction,'value') };

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%% Add axis checkbox
function CHECKaxis_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value of checkbox
handles.S.add_axis = get(handles.CHECKaxis,'value');

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%% Add colorbar checkbox
function CHECKcolorbar_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value of checkbox
handles.S.add_colorbar = get(handles.CHECKcolorbar,'value');

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%% Camlight pushbutton
function PUSHcamlight_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Modify light position vector
handles.S.light(handles.S.current_light).light_position = campos;
set(handles.H.light_handle(handles.S.current_light),...
    'position',handles.S.light(handles.S.current_light).light_position);

%Update GUI
guidata(hObject,handles);
spherium_update_gui_from_S

%%

%% Load settings pushbutton
function PUSHloadsettings_Callback(hObject, eventdata, handles)

%Return path to valid .mat file obtained via Explorer navigation by user
%which contains spherium setup information.
handles = guidata(hObject);
[filename, pathname, filterindex] = uigetfile([handles.last_lookin_dir,'\*.mat'], 'Choose a .mat file containing spherium settings');
if pathname ~=0
    handles.last_lookin_dir = pathname;
    
    %Load .mat file and update handles array of GUI with information contained
    %within the .mat file. If this doesn't work flag up an error dialogue box
    old_S = handles.S;
    try
        load( [pathname, '\', filename] );
        
        %Store spherium settings in handles array, to enable GUI callback
        %functions to see this data
        handles.S = S;
        
        %Update GUI
        guidata(hObject,handles);
        spherium_update_gui_from_S
        spherium_update;
    catch
        errordlg('.mat file selected does not contain the correct inputs for spherium','spherium error!')
        %Update GUI with original settings
        handles.S = old_S;
        guidata(hObject,handles);
        spherium_update_gui_from_S
        spherium_update;
        return
    end
end

%%

%% Load holes .jpg file %%
function PUSHloadholes_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
[filename, pathname, filterindex] = uigetfile([handles.last_lookin_dir,'\*.jpg'],...
    'Choose a .jpg file comprising the holes map');
if pathname ~=0
    handles.last_lookin_dir = pathname;
    copyfile( [pathname,'\',filename],'spherium\holes.jpg','f');
    guidata(hObject,handles);
    spherium_update;
end

%%

%% Load texture .jpg file %%
function PUSHloadtexture_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
[filename, pathname, filterindex] = uigetfile([handles.last_lookin_dir,'\*.jpg'],...
    'Choose a .jpg file comprising the texture map');
if pathname ~=0
    handles.last_lookin_dir = pathname;
    copyfile( [pathname,'\',filename],'spherium\texture.jpg','f');
    guidata(hObject,handles);
    spherium_update;
end

%%

% Load transparency .jpg file
function PUSHloadtransparency_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
[filename, pathname, filterindex] = uigetfile([handles.last_lookin_dir,'\*.jpg'],...
    'Choose a .jpg file comprising the transparency map');
if pathname ~=0
    handles.last_lookin_dir = pathname;
    copyfile( [pathname,'\',filename],'spherium\transparency.jpg','f');
    guidata(hObject,handles);
    spherium_update;
end

%%

%% Save settings pushbutton
function PUSHsavesettings_Callback(hObject, eventdata, handles)

%Return path to valid .mat file obtained via Explorer navigation by user
%which will contains spherium setup information.
[filename, pathname, filterindex] = uiputfile([handles.last_lookin_dir,'\*.mat'],...
    'Save spherium GUI settings in a .mat file');

%Save .mat file containing GUI settings
if ~isequal(filename,0) && ~isequal(pathname,0)
    handles = guidata(hObject);
    S = handles.S;
    save( [pathname,'\',filename], 'S' );
end

%%

%% Zoom in pushbutton +
function PUSHzoomin_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camzoom(1.2);
spherium_fast_update;

%%

%% Zoom out pushbutton -
function PUSHzoomout_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camzoom(1/1.2);
spherium_fast_update;

%%

%% Dolly up pushbutton
function PUSHup_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camdolly(0,-0.02,0);
spherium_fast_update;

%%

%% Dolly down pushbutton
function PUSHdown_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camdolly(0,0.02,0);
spherium_fast_update;

%%

%% Dolly right pushbutton
function PUSHright_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camdolly(-0.02,0,0);
spherium_fast_update;

%%

%% Dolly left pushbutton
function PUSHleft_Callback(hObject, eventdata, handles)
axes(handles.AXESspherium);
camdolly(0.02,0,0);
spherium_fast_update;

%%

%% Toggle camera roll %%
function TOGGLEroll_Callback(hObject, eventdata, handles)

%Make sure camera pan is untoggled
set( handles.TOGGLEpan, 'value', 0 );

%%

%% Toggle camera pan %%
function TOGGLEpan_Callback(hObject, eventdata, handles)

%Make sure camera roll is untoggled
set( handles.TOGGLEroll, 'value', 0 );

%%

%% Window Button Up function - complete change of light source angle with mouse movement
function windowbuttonup(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Set mouse button click status
handles.button='not clicked';

%Update GUI
guidata(gcf,handles);
spherium_fast_update;

%%

%% Window Button Down function - toggle change of light source angle with mouse movement
function windowbuttondown(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Set flag in handles structure which defines the type of click event.
handles.button = get(gcf, 'SelectionType');

%Double left mouse button click
if strcmp(handles.button,'open')
    handles.button = 'double click';
    
    %Single left mouse click
elseif strcmp(handles.button,'normal')
    handles.button = 'single click';
    
    %Shift + click or right & left mouse click
elseif strcmp(handles.button,'extend')
    handles.button = 'shift & click';
    
    %Control + click or right mouse click
elseif strcmp(handles.button,'alt')
    handles.button = 'ctrl & click';
else
    handles.button='not clicked';
end

%Store 'click point' coordinates (pixels from bottom left of the figure)
cp=get(handles.FIGUREspherium,'currentpoint');
handles.click_point_x = cp(1,1);
handles.click_point_y = cp(1,2);

%Get light azimuth and elevation when mouse was depressed (!). Note azimuth
%is measured anticlockwise from the -y axis (!!)
[theta,elev,range] = cart2sph( handles.S.light(handles.S.current_light).light_position(1),...
    handles.S.light(handles.S.current_light).light_position(2),...
    handles.S.light(handles.S.current_light).light_position(3) );
handles.light_azi_rad = theta + pi/2;
handles.light_elev_rad = elev;

%Update guidata
guidata(gcf,handles);

%%

%% Window Button Motion function - use mouse to define light source angle
%  or camera orbit for 3D rotation
function windowbuttonmotion(hObject, eventdata, handles)

if isfield(handles,'GUIisopen');
    if (gcf==handles.FIGUREspherium) && ( ~strcmp( handles.button,'not clicked') )
        
        %Get GUI data
        handles = guidata(gcf);
        
        %Store 'click point' coordinates (pixels from bottom left of the figure)
        cp=get(handles.FIGUREspherium,'currentpoint');
        x = cp(1,1);
        y = cp(1,2);
        
        %Get size of figure in pixels
        p = get( handles.FIGUREspherium,'position' );
        
        %Determine azimuth and elevation change from mouse movement
        u = 2 * (x - handles.click_point_x )/p(3) ;
        if u < -1
            u = -1;
        elseif u > 1
            u = 1;
        end
        az = 2*asin(u);
        v = 2 * (y - handles.click_point_y )/p(4) ;
        if v < -1
            v = -1;
        elseif u > 1
            v = 1;
        end
        el = 2*asin(v);
        
        %Only modify light if light toggle button is on
        if handles.S.light_toggle == 1
            %Modify light source angle depending of location of mouse button
            %within axes
            az = az + handles.light_azi_rad;
            el = el + handles.light_elev_rad;
            handles.S.light(handles.S.current_light).light_position =...
                norm( handles.S.light(handles.S.current_light).light_position ) *...
                [cos(el)*sin(az),-cos(el)*cos(az),sin(el) ];
            
            %Modify lighting
            axes( handles.AXESspherium );
            set(handles.H.light_handle(handles.S.current_light),...
                'position',handles.S.light(handles.S.current_light).light_position );
            
            %Update guidata
            guidata(gcf,handles);
        else
            %Otherwise rotate the axes via changing the camera viewpoint
            if get( handles.TOGGLEpan,'value')==0 && get( handles.TOGGLEroll,'value')==0
                
                %Orbit the camera (default)
                if strcmp(get(gcf, 'SelectionType'),'normal')
                    camorbit( -5*u, -5*v, 'camera' );
                else
                    %Holding down extra buttons rotates about the data
                    camorbit( -5*u, -5*v, 'data' );
                end
            elseif get( handles.TOGGLEpan,'value')==1
                
                %Pan the camera
                campan( -u/5, -v/5, 'camera' );
            elseif get( handles.TOGGLEroll,'value')==1
                
                %Roll the camera
                if v*u>0
                    camroll( (3*v+u) );
                else
                    camroll( (3*v+u) );
                end
            end
        end
    end
end

%%

%% Toggle auto axes scaling %%
function TOGGLEautoaxes_Callback(hObject, eventdata, handles)
spherium_update;

%%

%% Surface material popupmenu
function POPUPMENUmaterial_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Get string from menu
v = get( handles.POPUPMENUmaterial, 'value' );
handles.S.material = handles.H.materials{v};
material(handles.S.material);

%Update guidata
guidata(gcf,handles);

%%

%% Light color pushbutton
function PUSHlightcolor_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Determine light color from color picker
handles.S.light(handles.S.current_light).light_color =...
    uisetcolor( handles.S.light(handles.S.current_light).light_color, 'Set light colour' );
set( handles.PUSHlightcolor,'backgroundcolor',handles.S.light(handles.S.current_light).light_color );
set( handles.H.light_handle(handles.S.current_light),...
    'color',handles.S.light(handles.S.current_light).light_color );

%Update guidata
guidata(gcf,handles);

%%

%% Lighting style popupmenu
function POPUPMENUlightingstyle_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Get lighting style type from menu and update properties of Light # v
v = get( handles.POPUPMENUlightingstyle, 'value' );
handles.S.light(handles.S.current_light).light_style = handles.H.light_styles{v};
if strcmp( handles.S.light_styles{v}, 'none' )
    set( handles.H.light_handle(handles.S.current_light),'visible','off' );
else
    set( handles.H.light_handle(handles.S.current_light),'style',...
        handles.S.light(handles.S.current_light).light_style,'visible','on' );
end

%Update guidata
guidata(gcf,handles);

%%

%% Light model popupmenu
function POPUPMENUlightingmodel_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Get string from menu
v = get( handles.POPUPMENUlightingmodel, 'value' );
handles.S.lighting_model = handles.H.lighting_models{v};
lighting(handles.S.lighting_model);

%Update guidata
guidata(gcf,handles);

%%

%% Light toggle button.
function TOGGLElight_Callback(hObject, eventdata, handles)

%Get GUI data
handles = guidata(gcf);

%Update light toggle status
handles.S.light_toggle = get( handles.TOGGLElight,'value');

%Modify 3D rotation mode depending on toggle
axes( handles.AXESspherium );
if handles.S.light_toggle == 1
    handles.S.status = spherium_status('light toggled');
else
    handles.S.status = spherium_status('default');
end
set( handles.TEXTstatus,'string',handles.S.status );

%Update guidata
guidata(gcf,handles);

%%

%Load transparency map
function CHECKtransparency_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value of checkbox
handles.S.transparency = get(handles.CHECKtransparency,'value');

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%Load texture map
function CHECKtexture_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value of checkbox
handles.S.texture = get(handles.CHECKtexture,'value');

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%Load holes map
function CHECKholes_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Get value of checkbox
handles.S.holes = get(handles.CHECKholes,'value');

%Update GUI
guidata(hObject,handles);
spherium_update;

%%

%Select light popupmenu
function POPUPMENUselectlight_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Modify current light selected
handles.S.current_light = get(handles.POPUPMENUselectlight,'value' );

%Update GUI
guidata(hObject,handles);
spherium_update_gui_from_S

%%

%Light x coordinate edit box
function EDITlightx_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Update handles array based upon new light position coordinates
handles.S.light( handles.S.current_light ).light_position(1) = ...
    str2num( get( handles.EDITlightx,'string' ) );

%Modify lighting
axes( handles.AXESspherium );
set(handles.H.light_handle(handles.S.current_light),...
    'position',handles.S.light(handles.S.current_light).light_position );

%Update GUI
guidata(hObject,handles);

%%

%Light y coordinate edit box
function EDITlighty_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Update handles array based upon new light position coordinates
handles.S.light( handles.S.current_light ).light_position(2) = ...
    str2num( get( handles.EDITlighty,'string' ) );

%Modify lighting
axes( handles.AXESspherium );
set(handles.H.light_handle(handles.S.current_light),...
    'position',handles.S.light(handles.S.current_light).light_position );

%Update GUI
guidata(hObject,handles);

%%

%Light z coordinate edit box
function EDITlightz_Callback(hObject, eventdata, handles)

%Get data from GUI
handles = guidata(hObject);

%Update handles array based upon new light position coordinates
handles.S.light( handles.S.current_light ).light_position(3) = ...
    str2num( get( handles.EDITlightz,'string' ) );

%Modify lighting
axes( handles.AXESspherium );
set(handles.H.light_handle(handles.S.current_light),...
    'position',handles.S.light(handles.S.current_light).light_position );

%Update GUI
guidata(hObject,handles);

%%

%% Close function
function spherium_close(hObject, eventdata, handles)

%Saves settings (in structure handles.S) to a .mat file
handles = guidata(hObject);
S = handles.S;
try
    save spherium\last_spherium_settings S
    last_lookin_dir = handles.last_lookin_dir;
    save spherium\spherium_store last_lookin_dir
catch err
    %Possible error is one may not have write-access to the current drive
    errordlg( {'Have you got write access to the current drive?',
        'If not spherium will not be able to save your settings on exiting.'},...
        'spherium error' );
end

%Close previous GUI
delete( findobj('tag',['FIGURE',handles.S.spheria] ));

%Closes the spherium main GUI
delete( handles.FIGUREspherium );

%Remove current directory tree from MATLAB path
rmpath(genpath(pwd))

%End of code







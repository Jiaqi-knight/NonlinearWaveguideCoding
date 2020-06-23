%spherium_defaults
%
% S  Structure with parameters that the user will alter
% H  Data structure internal to spherium. i.e. user cannot change via
%    loading a file.
function [S,H] = spherium_defaults

%Plot view azimuth and elevation /deg
S.camera_position = [15.6297 6.2151 -4.4025];
S.camera_target = [0.103866 -0.104356 0.180768];
S.camera_view_angle = 4.83194;
S.camera_up_vector = [0.143934 0.324424 0.934897];
[az,el,R] = cart2sph( S.camera_position(1),...
    S.camera_position(2),...
    S.camera_position(3) );
S.azi_deg = az*180/pi;
S.elev_deg = el*180/pi;

%Output image size and resolution in dots per inch
S.output_size = 'A4 landscape';
H.output_sizes = {
    'A6 landscape',...
    'A6 portrait',...
    'A5 landscape',...
    'A5 portrait',...
    'A4 landscape',...
    'A4 portrait',...
    'A3 landscape',...
    'A3 portrait',...
    'A2 landscape',...
    'A2 portrait',...
    'A1 landscape',...
    'A1 portrait',...
    'A0 landscape',...
    'A0 portrait',...
    '6" x 4"',...
    '7" x 5"',...
    '8" x 6"',...
    '10" x 12"',...
    '10" x 15"'
    };
S.output_DPI = 300;

%Spheria and their default parameters
S.spheria = 'ammonite';
[S.spherias, S.SPHERIA] = spherium_get_spheria_defaults( 'spherium\SPHERIA' );

%Plot rendering. Note OpenGL is usually the best but may not be avaliable
%on all platforms
S.renderer = 'OpenGL';
H.renderers = {
    'OpenGL',...
    'zbuffer',...
    };

%Surface material
S.material = 'default';
H.materials = {
    'default',...
    'shiny',...
    'dull',...
    'metal'
    };

%Transparency?
S.transparency = 0;

%Holes?
S.holes = 0;

%Texture map?
S.texture = 0;

%Number of surface points per dimension (total is S.N ^ 2)
S.N = 500;

%Add x,y,z axis to plot
S.add_axis = 0;

%Add colorbar
S.add_colorbar = 0;
H.colorbar_handle = [];

%Colour maps
H.colour_maps = {
    'hsv',...       %Hue-saturation-value color map.
    'hot',...       %Black-red-yellow-white color map.
    'gray',...      %Linear gray-scale color map.
    'bone',...      %Gray-scale with tinge of blue color map.
    'copper',...    %Linear copper-tone color map.
    'pink',...      %Pastel shades of pink color map.
    'white',...     %All white color map. Use with a non-white NaN colour
    'flag',...      %Alternating red, white, blue, and black color map.
    'lines',...     %Color map with the line colors.
    'colorcube',... %Enhanced color-cube color map.
    'vga',...       %Windows colormap for 16 colors.
    'jet',...       %Variant of HSV.
    'black & white',... %Alternating black and white - similar to 'prism'                   
    'prism',...     %Prism color map.
    'cool',...      %Shades of cyan and magenta color map.
    'autumn',...    %Shades of red and yellow color map.
    'spring',...    %Shades of magenta and yellow color map.
    'winter',...    %Shades of blue and green color map.
    'summer',...    %Shades of green and yellow color map.
    'bespoke'};

S.colour_map = 'jet';

%Colour functions
S.colorfunction = 'None';
H.colorfunctions = {
    'None',...
    'Log',...
    'Sine',...
    'Sine 10',...
    'Sine 200'
    };

%Default output filename
S.filename = 'spherium image.png';

%Default status box message
S.status = spherium_status('default');

%Lighting settings ....

%Lighting algorithm
S.lighting_model = 'phong';
H.lighting_models = {
    'phong',...
    'gouraud',...
    'flat',...
    'none'
    };

%Light styles
H.light_styles = {
    'infinite',...
    'local',...
    'none'
    };

%Light 3D rotation toggle flag
S.light_toggle = 0;

%List of lights
H.lights = {
    'Light 1',...
    'Light 2',...
    'Light 3'
    };
S.current_light = 1;

%Lighting color
S.light(1).light_color = [1 1 1];
S.light(2).light_color = [1 1 1];
S.light(3).light_color = [1 1 1];

%Light style
S.light(1).light_style = 'local';
S.light(2).light_style = 'none';
S.light(3).light_style = 'none';

%Light position 
S.light(1).light_position = 1e3*[ 1.0080, 0.0875, 0.0258];
S.light(2).light_position = 1e3*[ 0, 0, 0];
S.light(3).light_position = 1e3*[ -1, 0.1, -0.3];

%

%Handle for surface plot
H.surface_handle = [];

%Handle for spiral plot figure
H.ammonite_spiral_plot = [];

%End of code
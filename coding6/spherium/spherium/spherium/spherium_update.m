%spherium_update
% Updates the spherium surface generic features such as lighting, holes,
% transparency and textures.

function spherium_update

%Get handle to main GUI figure and then extract gui data structure S
mainfig = findobj('tag','FIGUREspherium');
h = guidata( mainfig );

%Point to main plot axes
figure(h.FIGUREspherium);
axes(h.AXESspherium);

%Delete any previous surfaces or colorbars
delete( h.H.colorbar_handle ); h.H.colorbar_handle = [];
delete( h.H.surface_handle ); h.H.surface_handle = [];
cla;
guidata( mainfig, h );

%Re-draw spheria from '...._update_surface function associated with
%spheria'
eval([h.S.spheria,'_update_surface;']);

%Set generic axes properties based upon parameters stored in the GUI data
spherium_update_gui_from_S;
h = guidata( mainfig );

%

% NORMALIZE
x = get( h.H.surface_handle, 'Xdata');
y = get( h.H.surface_handle, 'Ydata');
z = get( h.H.surface_handle, 'Zdata');
r = sqrt( x.^2 + y.^2 + z.^2 );
norm = max(max(r));
r = r/norm;
x = x/norm;
y = y/norm;
z = z/norm;

%

% COLOUR FUNCTIONS

%Get colour data matrix and normalize
C = get( h.H.surface_handle, 'CData' );
C = C/max(max(C));

% Define colour function operator
if strcmp( h.S.colorfunction, 'Sine 200' )
    colour_func = 'sin(200*pi*C)';
elseif strcmp( h.S.colorfunction, 'Sine 10' )
    colour_func = 'sin(20*pi*C)';
elseif strcmp( h.S.colorfunction, 'Sine' )
    colour_func = 'sin(pi*C)';
elseif strcmp( h.S.colorfunction, 'Log' )
    colour_func = 'log(0.1 + abs(C))';
else
    colour_func='C';
end
eval( ['C=',colour_func,';'] );
if max(max(C)) - min(min(C)) > 0
    caxis( [min(min(C)),max(max(C))] );
end
set( h.H.surface_handle, 'CData', C );

%

% COLORMAP

%Make new interpolated colormap
N = round(sqrt(numel(C)));
if strcmp(h.S.colour_map,'black & white')==1
    %Create alternating black & white colormap
    map=zeros(N,3);
    for k=1:N
        if mod(k,2)==0
            map(k,:)=[1,1,1];
        end
    end
elseif strcmp(h.S.colour_map,'bespoke')==1
    %Get current colormap
    map = colormap;
else
    map = colormap(h.S.colour_map);
end
colormap(map);
interp_colormap(N);

%

% COLORBAR
if h.S.add_colorbar==1
    h.H.colorbar_handle = colorbar('fontsize',8,...
        'location','North',...
        'position',[0.02 0.95 0.4 0.03]);
end

%

% ADD AXIS %

%Add x,y,z lines
hold on
if h.S.add_axis==1
    line([0,1.5*max(max(abs(x)))],[0,0],[0,0],'color','r');
    line([0,0],[0,1.5*max(max(abs(y)))],[0,0],'color','g');
    line([0,0],[0,0],[0,1.5*max(max(abs(z)))],'color','b');
end

%

% HOLES MAP

%If holes map option is chosen, punch holes in surface.
if h.S.holes == 1
    
    %Load image of texture from file texture.jpg
    [I,map] = imread( 'spherium\holes.jpg' );
    
    %Get size of surface matrix
    C = get( h.H.surface_handle, 'CData');
    [M,N,P] = size(C);
    
    %Make sure it is only one 'page' i.e. greyscale. If not mean over the
    %other colours.
    [r,c,p] = size(I);
    if p~=1
        II = zeros(r,c);
        II = (1/3) * ( I(:,:,1)+I(:,:,2)+I(:,:,3) );
        I = II;
    end
    clear II
    
    %Interpolate image
    I = interp_image( I , N , M , '' ) ;
    
    %Convert into range 0-1 and flip
    I=I-min(min(I));
    I = I/max(max(I));
    I = flipdim(I,1);
    
    %Convert to Black and White
    I(I>=0.5) = 1;
    I(I<0.5) = 0;
    
    %Set NaN values at hole positions (black in the hole map image)
    x(I==0) = NaN;
    y(I==0) = NaN;
    z(I==0) = NaN;
    set(  h.H.surface_handle, 'Xdata',x );
    set(  h.H.surface_handle, 'Ydata',y );
    set(  h.H.surface_handle, 'Zdata',z );
    clear X Y Z
end

%

% TEXTURE MAP

%If texture map option is chosen, change colouring of surface.
if h.S.texture == 1
    
    %Load image of texture from file texture.jpg
    [I,map] = imread( 'spherium\texture.jpg' );
    
    %Get size of surface matrix
    C = get( h.H.surface_handle, 'CData');
    [M,N,P] = size(C);
    
    %Interpolate image
    I = interp_image( I , N , M , '' ) ;
    
    %Convert into range 0-1 and flip
    I=I/255;
    I = flipdim(I,1);
    
    %texture map surface
    set( h.H.surface_handle ,'CDatamapping', 'Direct','CData', I);
    
end

%

% TRANSPARENCY MAP

%If transparency option is chosen, change colouring of surface.
if h.S.transparency ==1
    
    %Load image of texture from file texture.jpg
    [I,map] = imread( 'spherium\transparency.jpg' );
    
    %Get size of surface matrix
    C = get( h.H.surface_handle, 'CData');
    [M,N,P] = size(C);
    
    %Interpolate image
    I = interp_image( I , N , M , '' ) ;
    
    %Make sure it is only one 'page' i.e. greyscale. If not mean over the
    %other colours.
    [r,c,p] = size(I);
    if p~=1
        II = zeros(r,c);
        II = (1/3) * ( I(:,:,1)+I(:,:,2)+I(:,:,3) );
        I = II;
    end
    clear II
    
    %Convert into range 0-1 and flip
    I=I-min(min(I));
    I = I/max(max(I));
    I = flipdim(I,1);
    
    %Modify alpha map of surface
    set( h.H.surface_handle ,'AlphaDataMapping','none','AlphaData', I,...
        'edgealpha','interp','facealpha','interp');
end

%

%LIGHTING

%Create light objects
for n=1:length(h.S.lights)
    if strcmp( h.S.light(n).light_style,'none' )
        h.H.light_handle(n) = light( 'position',h.S.light(n).light_position,...
            'color',h.S.light(n).light_color,'style','local','visible','off' );
    else
        h.H.light_handle(n) = light( 'position',h.S.light(n).light_position,...
            'color',h.S.light(n).light_color,'style',h.S.light(n).light_style );
    end
end

%Set surface material lighting
material( h.S.material );
lighting( h.S.lighting_model );

%

%Update gui data
guidata( mainfig, h );

%End of code

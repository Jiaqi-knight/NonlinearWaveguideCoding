%spherium_update_gui_from_S
% Updates the spherium GUI basd upon settings stored within structure S
% held within the GUI data structure.

function spherium_update_gui_from_S

%Get handle to main spherium GUI figure
mainfig = findobj('tag','FIGUREspherium');
h = guidata( mainfig );

%Update status box
set( h.TEXTstatus,'string',h.S.status );

%Update version
set( h.TEXTversion,'string',spherium_version );

%Axes properties: Get these from the current spheria. During startup with
%default parameters there will be no axes data, hence the empty check. The
%default axes_properties arrays will be [].
eval( ['emptycheck = isempty(h.S.SPHERIA.',h.S.spheria,'.axes_properties);']);
if emptycheck == 0
    eval( ['set( h.AXESspherium,h.S.SPHERIA.',h.S.spheria,...
        '.axes_properties_fields,h.S.SPHERIA.',h.S.spheria,'.axes_properties );'] )
end

%Plot view azimuth and elevation /deg
set( h.EDITazi,'string', num2str(h.S.azi_deg) ) ;
set( h.EDITelev,'string', num2str(h.S.elev_deg) );
set( h.SLIDERazi,'value', (h.S.azi_deg+180)/360 );
set( h.SLIDERelev,'value', (h.S.elev_deg+90)/180 );

%Output image size and resolution in dots per inch
set( h.POPUPMENUoutputsize,'string',h.H.output_sizes );
set( h.POPUPMENUoutputsize,'value', strmatch( h.S.output_size, h.H.output_sizes, 'exact' ));
set( h.EDITDPI,'string',num2str(h.S.output_DPI) );

%Plot type
set( h.POPUPMENUspheria, 'string', h.S.spherias );
set( h.POPUPMENUspheria,'value',...
    strmatch( h.S.spheria, h.S.spherias, 'exact' ));

%Plot rendering. Note OpenGL is usually the best but may not be avaliable
%on all platforms
set( h.POPUPMENUrenderer, 'string', h.H.renderers );
set( h.POPUPMENUrenderer,'value',...
    strmatch( h.S.renderer, h.H.renderers, 'exact' ));

%Surface plot material
set( h.POPUPMENUmaterial, 'string', h.H.materials );
set( h.POPUPMENUmaterial,'value',...
    strmatch( h.S.material, h.H.materials, 'exact' ));

%Transparency
set( h.CHECKtransparency,'value',h.S.transparency );

%Holes
set( h.CHECKholes,'value',h.S.holes );

%Add x,y,z axis to plot
set( h.CHECKaxis,'value',h.S.add_axis );

%Add colorbar
set( h.CHECKcolorbar,'value',h.S.add_colorbar );

%General colouring settings ...
set( h.POPUPMENUcolormap,'string',h.H.colour_maps );
set( h.POPUPMENUcolormap,'value',...
    strmatch( h.S.colour_map , h.H.colour_maps, 'exact' ));
set( h.POPUPMENUcolorfunction,'string',h.H.colorfunctions );
set( h.POPUPMENUcolorfunction,'value',...
    strmatch( h.S.colorfunction , h.H.colorfunctions, 'exact' ));

%Lighting settings ....

%List of lights
set (h.POPUPMENUselectlight,'string',h.H.lights );
set( h.POPUPMENUselectlight,'value',h.S.current_light);

%Lightning model
set( h.POPUPMENUlightingmodel, 'string', h.H.lighting_models );
set( h.POPUPMENUlightingmodel,'value',...
    strmatch( h.S.lighting_model, h.H.lighting_models, 'exact' ));

%Lighting style
set( h.POPUPMENUlightingstyle, 'string', h.H.light_styles );
set( h.POPUPMENUlightingstyle,'value',...
    strmatch( h.S.light(h.S.current_light).light_style, h.H.light_styles, 'exact' ));

%Color of light color pushbutton
set( h.PUSHlightcolor,'backgroundcolor',h.S.light(h.S.current_light).light_color );

%Light position edit boxes
set( h.EDITlightx,'string',num2str( h.S.light(h.S.current_light).light_position(1) ) );
set( h.EDITlighty,'string',num2str( h.S.light(h.S.current_light).light_position(2) ) );
set( h.EDITlightz,'string',num2str( h.S.light(h.S.current_light).light_position(3) ) );

%Light toggle status
set( h.TOGGLElight,'value',h.S.light_toggle );

%Update GUI
guidata( mainfig, h );

%End of code
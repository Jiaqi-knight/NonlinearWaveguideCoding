%spherium_fast_update

function spherium_fast_update

%Get GUI data
h = guidata(gcf);

%Determines current azimuth and elevation of main Spherium axes and updates
%sliders and edit boxes. It is assumed these have been modified via 3D
%rotation.
p = get(h.AXESspherium,'CameraPosition');
[az,el,R] = cart2sph( p(1),p(2),p(3) );
h.S.azi_deg = az*180/pi;
h.S.elev_deg = el*180/pi;
S.camera_position = p;
S.camera_target = get(h.AXESspherium,'CameraPosition');
S.camera_view_angle = get(h.AXESspherium,'CameraViewAngle');
S.camera_up_vector = get(h.AXESspherium,'CameraUpVector');

%Store current axes properties in the GUI data of the current SPHERIA. This
%means each spheria can have there own display properties rather than adopt
%the properties of the previous spheria.
eval(['[h.S.SPHERIA.',h.S.spheria,'.axes_properties_fields,h.S.SPHERIA.',...
    h.S.spheria,'.axes_properties] = getax(h.AXESspherium);']);

%Update GUI
guidata( gcf, h);
spherium_update_gui_from_S;

%End of code
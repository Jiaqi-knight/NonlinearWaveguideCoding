%spherium_status
% Function which returns a text string corresponding to the status of spherium.m

function s = spherium_status( status )

if strcmp( status, 'default' )
    s = ['Welcome to Spherium. Dragging the mouse in the main axes will result in a 3D rotation.',...
    'Use the + and - buttons to zoom in and out, and the >, < etc to translate the figure. ',...
    'Hydrogenic orbital spheria take the form H XN e.g. H P1. Note the blue square must be pressed ',...
    'to update Spherium following 3D rotation.'];
elseif strcmp( status, 'light toggled' )
    s = ['Light mode toggled. Dragging the mouse in the main axes will result in a .',...
    '3D rotation of the light source. The + and - buttons can be used to bring the source ',...
    'closer or further away from the spheria.'];
elseif strcmp( status, 'Saving PNG' )
    s = ['Saving PNG image of spheria. A new figure window may pop up during this operation. ',...
        'Please wait until the image creation process has completed.' ] ;
end

%End of code
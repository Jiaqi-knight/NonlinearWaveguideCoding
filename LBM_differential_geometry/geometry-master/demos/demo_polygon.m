%DEMO_POLYGON  One-line description here, please.
%
%   output = demo_polygon(input)
%
%   Example
%   demo_polygon
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-29,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.


%% Create polygon

% create a polygon with arbitrary vertices
verts = [...
    10 60; 20 40; 30 40; 20 10; 70 40; 80 20; ...
    90 50; 60 80; 50 60; 40 80; 30 60];
poly = Polygon2D(verts);

% display polygon
figure; hold on; 
draw(poly);
axis equal; axis([0 100 0 100]);
drawnow;


%% Create other geometries from polygon

% draw bounding box
bbox = boundingBox(poly);
draw(bbox, 'k')

% draw vertices
verts = vertices(poly);
draw(verts, 'ks');

% draw the centroid
cent = centroid(poly);
draw(cent, 'bo');


%% Compute area and perimeter

a = area(poly);
p = perimeter(poly);

fprintf('Area = %7.2f\n', a);
fprintf('Perimeter = %7.2f\n', p);

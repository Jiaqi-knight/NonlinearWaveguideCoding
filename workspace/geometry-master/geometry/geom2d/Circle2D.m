classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Circle2D < Curve2D
% A circle in the plane.
%
%   Class Circle2D
%
%   Example
%     Circle2D
%
%   See also
%     Ellipse2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-24,    using Matlab 9.6.0.1072779 (R2019a)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    CenterX = 0;
    CenterY = 0;
    Radius = 1;
    
end % end properties


%% Constructor
methods
    function obj = Circle2D(varargin)
    % Constructor for Circle2D class.

        switch nargin
            case 0
                % nothing to do
            case 1
                var1 = varargin{1};
                if size(var1, 2) ~= 3
                    error('Creating a circle requires an array with three columns, not %d', size(var1, 2));
                end
                obj.CenterX = var1(1);
                obj.CenterY = var1(2);
                obj.Radius = var1(3);
            case 2
                var1 = varargin{1};
                if size(var1, 2) ~= 2
                    error('Creating a circle requires an array with two columns, not %d', size(var1, 2));
                end
                obj.CenterX = var1(1);
                obj.CenterY = var1(2);
                obj.Radius = varargin{2};
        end
    end

end % end constructors


%% Methods specific to Circle2D
methods
    function center = center(obj)
        % Returns the center of this circle as a Point2D.
        center = Point2D(obj.CenterX, obj.CenterY);
    end
    
    function poly = asPolyline(obj, varargin)
        % Convert this circle into a polyline.
        
        % determines number of points
        N = 64;
        if ~isempty(varargin)
            N = varargin{1};
        end
        
        % create circle
        t = linspace(0, 2*pi, N+1)';
        t(end) = [];
        
        % coordinates of circle points
        x = obj.CenterX + obj.Radius * cos(t);
        y = obj.CenterY + obj.Radius * sin(t);

        poly = LinearRing2D([x y]);
    end
end

%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transform) %#ok<STOUT>
        % Apply a geometric transform to this geometry.
        error('Transform not implemented for Circles');
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this geometry.
        extX = [obj.CenterX - obj.Radius obj.CenterX + obj.Radius];
        extY = [obj.CenterY - obj.Radius obj.CenterY + obj.Radius];
        box = Box2D([extX extY]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % compute a set of coordinates for drawing circle
        N = 72;
        t = linspace(0, 2*pi, N+1);
        xt = obj.CenterX + obj.Radius * cos(t);
        yt = obj.CenterY + obj.Radius * sin(t);
        
        hh = plot(ax, xt, yt, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        res = Circle2D([obj.CenterX obj.CenterY obj.Radius] * factor);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        res = Circle2D([obj.CenterX obj.CenterY] + shift,  obj.Radius);
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this circle.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.
        
        center2 = rotate(center(obj), angle, varargin{:});
        res = Circle2D([center2.X center2.Y obj.Radius]);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Circle2D', 'CenterX', obj.CenterX, 'CenterY', obj.CenterY, 'Radius', obj.Radius);
    end
end
methods (Static)
    function circ = fromStruct(str)
        % Create a new instance from a structure.
        circ = Circle2D([str.CenterX str.CenterY str.Radius]);
    end
end

end % end classdef


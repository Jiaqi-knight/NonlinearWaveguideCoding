classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LineSegment2D < Curve2D
% A line segment defined by its two extremities.
%
%   Class LineSegment2D
%
%   Example
%     P1 = Point2D(20, 10);
%     P2 = Point2D(40, 20);
%     L = LineSegment2D(P1, P2);
%     draw(L);
%
%   See also
%     Point2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The x-coordinate of the source point.
    X1 = 0;
    % The y-coordinate of the source point.
    Y1 = 0;

    % The x-coordinate of the target point.
    X2 = 1;
    % The y-coordinate of the target point.
    Y2 = 0;
    
end % end properties


%% Constructor
methods
    function obj = LineSegment2D(varargin)
        % Constructor for LineSegment2D class
        
        if nargin == 0
            % Default constructor: unit line segment
            
        elseif nargin == 1
            % Copy constructor
            if ~isa(varargin{1}, 'LineSegment2D')
                error('Requires a LineSegment2D as input');
            end
            var1 = varargin{1};
            obj.X1 = var1.X1;
            obj.Y1 = var1.Y1;
            obj.X2 = var1.X2;
            obj.Y2 = var1.Y2;
            
        elseif nargin == 2
            p1 = varargin{1};
            if isa(p1, 'Point2D')
                obj.X1 = p1.X;
                obj.Y1 = p1.Y;
            else
                obj.X1 = p1(1);
                obj.Y1 = p1(2);
            end
              
            p2 = varargin{2};
            if isa(p2, 'Point2D')
                obj.X2 = p2.X;
                obj.Y2 = p2.Y;
            else
                obj.X2 = p2(1);
                obj.Y2 = p2(2);
            end
        end

    end

end % end constructors


%% Methods specific to LineSegment2D
methods
    function [dist, pos] = distancePoint(obj, point)
        % Minimum distance between point(s) and this line segment.
        
        % direction vector of the line segment
        vx = obj.X2 - obj.X1;
        vy = obj.Y2 - obj.Y1;
        
        % squared length of edges, with a check of validity
        delta = vx .* vx + vy .* vy;
        
        % difference of coordinates between points and first point
        dx  = point(:, 1) - obj.X1;
        dy  = point(:, 2) - obj.Y1;
        
        % compute position of points projected on the supporting line,
        % by using normalised dot product (NP-by-NE array)
        if delta > eps
            % ensure projected point is located on the edge
            pos = min(max((dx * vx + dy * vy) / delta, 0), 1);
        else
            % consider point1 is the closest egde point
            pos = 0;
        end
        
        % compute distance between point and its projection on the edge
        dist = hypot(pos * vx - dx, pos * vy - dy);
    end
    
    function pm = middlePoint(obj)
        % Return the middle point of this line segment.
        xm = (obj.X1 + obj.X2) / 2;
        ym = (obj.Y1 + obj.Y2) / 2;
        pm = Point2D(xm, ym);
    end
end

%% Methods generic to curve objects
methods
    function l = length(obj)
        % Return the length of this line segment.
        dx = obj.X2 - obj.X1;
        dy = obj.Y2 - obj.Y1;
        l = sqrt(dx*dx + dy*dy);
    end
    
    function res = reverse(obj)
        % Reverse this line segment.
        res = LineSegment2D(obj.P2, obj.P1);
    end
    
    function p1 = firstPoint(obj)
        % Return the first point of this line segment.
        p1 = Point2D(obj.X1, obj.Y1);
    end
    
    function p2 = lastPoint(obj)
        % Return the last point of this line segment.
        p2 = Point2D(obj.X2, obj.Y2);
    end
end


%% Methods implementing the Geometry2D interface
methods
    function [dist, pos] = distance(obj, point)
        % Distance between point(s) and this line segment.
        %
        % Example:
        %   seg = LineSegment2D(Point2D(10, 10), Point2D(70, 55));
        %   p = Point2D(56, 32);
        %   distance(seg, p)
        %   ans =
        %       10
        %
        % see also
        %   distancePoint
        if isa(point, 'Point2D')
            [dist, pos] = distancePoint(obj, [point.X point.Y]);
        end
    end
    
    function res = transform(obj, transfo)
        % Apply a geometric transform to this line segment.
        p1t = transformPoint(transfo, [obj.X1 obj.Y1]);
        p2t = transformPoint(transfo, [obj.X2 obj.Y2]);
        res = LineSegment2D(p1t, p2t);
    end
    
    function box = boundingBox(obj)
        % Returns the bounding box of this geometry.
        x = sort([obj.X1 obj.X2]);
        y = sort([obj.Y1 obj.Y2]);
        box = Box2D([x y]);
    end
    
    function h = draw(varargin)
        % Draws the current geometry, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % plot line segment
        xdata = [obj.X1 obj.X2];
        ydata = [obj.Y1 obj.Y2];
        hh = plot(ax, xdata, ydata, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Returns a scaled version of this line segment.
        res = LineSegment2D([obj.X1 obj.Y1] * factor, [obj.X2 obj.Y2] * factor);
    end
    
    function res = translate(obj, shift)
        % Returns a translated version of this line segment.
        res = LineSegment3D([obj.X1 obj.Y1] + shift, [obj.X2 obj.Y2] + shift);
    end
    
    function res = rotate(obj, varargin)
        % Returns a rotated version of this line segment.
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        rot = createRotation(origin, deg2rad(angle));
        p1t = transformPoint(rot, [obj.X1 obj.Y1]);
        p2t = transformPoint(rot, [obj.X2 obj.Y2]);
        
        res =  LineSegment2D(p1t, p2t);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LineSegment2D', ...
            'X1', obj.X1, 'Y1', obj.Y1, ...
            'X2', obj.X2, 'Y2', obj.Y2);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        p1 = [str.X1 str.Y1];
        p2 = [str.X2 str.Y2];
        line = LineSegment2D(p1, p2);
    end
end

end % end classdef


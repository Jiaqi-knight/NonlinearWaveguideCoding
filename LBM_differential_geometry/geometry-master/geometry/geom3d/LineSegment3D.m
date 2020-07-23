classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LineSegment3D < Geometry3D
% A 3D line segment defined by its two extremities.
%
%   Class LineSegment3D
%
%   Example
%     P1 = Point3D(30, 20, 10);
%     P2 = Point3D(50, 40, 20);
%     L = LineSegment3D(P1, P2);
%     draw(L);
%
%   See also
%     Point3D, LineSegment2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-02-06,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The x-coordinate of the source point.
    X1 = 0;
    % The y-coordinate of the source point.
    Y1 = 0;
    % The z-coordinate of the source point.
    Z1 = 0;

    % The x-coordinate of the target point.
    X2 = 1;
    % The y-coordinate of the target point.
    Y2 = 0;
    % The z-coordinate of the target point.
    Z2 = 0;
    
end % end properties


%% Constructor
methods
    function obj = LineSegment3D(varargin)
        % Constructor for LineSegment3D class.
        
        if nargin == 0
            % Default constructor: unit line segment
            return;
            
        elseif nargin == 1
            % Copy constructor
            if ~isa(varargin{1}, 'LineSegment3D')
                error('Requires a LineSegment3D as input');
            end
            var1 = varargin{1};
            obj.X1 = var1.X1;
            obj.Y1 = var1.Y1;
            obj.Z1 = var1.Z1;
            obj.X2 = var1.X2;
            obj.Y2 = var1.Y2;
            obj.Z2 = var1.Z2;
            
        elseif nargin == 2
            p1 = varargin{1};
            if isa(p1, 'Point3D')
                obj.X1 = p1.X;
                obj.Y1 = p1.Y;
                obj.Z1 = p1.Z;
            else
                obj.X1 = p1(1);
                obj.Y1 = p1(2);
                obj.Z1 = p1(3);
            end
              
            p2 = varargin{2};
            if isa(p2, 'Point3D')
                obj.X2 = p2.X;
                obj.Y2 = p2.Y;
                obj.Z2 = p2.Z;
            else
                obj.X2 = p2(1);
                obj.Y2 = p2(2);
                obj.Z2 = p2(3);
            end
        end

    end

end % end constructors


%% Methods generic to curve objects
methods
    function l = length(obj)
        dx = obj.X2 - obj.X1;
        dy = obj.Y2 - obj.Y1;
        dz = obj.Z2 - obj.Z1;
        l = sqrt(dx*dx + dy*dy + dz*dz);
    end
    
    function res = reverse(obj)
        res = LineSegment3D([obj.X2 obj.Y2 obj.Z2], [obj.X1 obj.Y1 obj.Z1]);
    end
    
    function p1 = firstPoint(obj)
        p1 = Point3D(obj.X1, obj.Y1, obj.Z1);
    end
    
    function p2 = lastPoint(obj)
        p2 = Point3D(obj.X2, obj.Y2, obj.Z2);
    end
end


%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transfo)
        % Apply a geometric transform to this line segment.
        p1t = transformPoint(transfo, [obj.X1 obj.Y1 obj.Z1]);
        p2t = transformPoint(transfo, [obj.X2 obj.Y2 obj.Z2]);
        res = LineSegment3D(p1t, p2t);
    end
    
    function box = boundingBox(obj)
        % Returns the bounding box of this geometry.
        x = sort([obj.X1 obj.X2]);
        y = sort([obj.Y1 obj.Y2]);
        z = sort([obj.Z1 obj.Z2]);
        box = Box3D([x y z]);
    end
    
    function h = draw(varargin)
        % Draws the current geometry, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % plot line segment
        xdata = [obj.X1 obj.X2];
        ydata = [obj.Y1 obj.Y2];
        zdata = [obj.Z1 obj.Z2];
        hh = plot3(ax, xdata, ydata, zdata, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Returns a scaled version of this geometry.
        res = LineSegment3D([obj.X1 obj.Y1 obj.Z2] * factor, [obj.X2 obj.Y2 obj.ZZ] * factor);
    end
   
    function res = translate(obj, shift)
        % Returns a translated version of this geometry.       
        res = LineSegment3D([obj.X1 obj.Y1 obj.Z2] + shift, [obj.X2 obj.Y2 obj.ZZ] + shift);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LineSegment3D', ...
            'X1', obj.X1, 'Y1', obj.Y1, 'Z1', obj.Z1, ...
            'X2', obj.X2, 'Y2', obj.Y2, 'Z2', obj.Z2);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        p1 = [str.X1 str.Y1 str.Z1];
        p2 = [str.X2 str.Y2 str.Z2];
        line = LineSegment3D(p1, p2);
    end
end

end % end classdef


classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MultiCurve2D < Geometry2D
% A set of curves in the plane.
%
%   Class MultiCurve2D
%
%   Example
%     seg1 = LineSegment2D(Point2D(4, 5), Point2D(6,5));
%     seg2 = LineSegment2D(Point2D(5, 4), Point2D(5, 6));
%     circ = Circle2D([7, 3], 2);
%     mc = MultiCurve2D(seg1, seg2, circ);
%     figure; axis equal; axis([0 10 0 10]); hold on;
%     draw(mc, 'b');
%
%   See also
%     Curve2D, MultiPoint2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The inner array of curves, as a 1-by-N cell array containing Curve2D instances
    Curves;
    
end % end properties


%% Constructor
methods
    function obj = MultiCurve2D(varargin)
        % Constructor for MultiCurve2D class.
    
        if isempty(varargin)
            % Empty constuctor
            obj.Curves = {};
        elseif nargin == 1
            % either copy constructor or initialisation constructor
            var1 = varargin{1};
            if isa(var1, 'MultiCurve2D')
                % copy constructor
                % TODO: implement copy constructor 
                % (requires implementing Geometry.clone method)
                error('Not yet implemented...');
            
            elseif isa(var1, 'Curve2D')
                obj.Curves = {var1};
            else
                error('wrong input argument for creating a MultiCurve2D');
            end
        else
            % initialisation constructor
            obj.Curves = cell(1, nargin);
            for i = 1:nargin
                obj.Curves{i} = varargin{i};
            end
        end
            
    end

end % end constructors


%% Methods specific to MultiCurve2D
methods
end % end methods


%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transfo)
        % Applies a geometric transform to this geometry
        curves = cell(size(obj.Curves));
        for i = 1:numel(curves)
            curves{i} = transform(curves{i}, transfo);
        end
        res = MultiCurve2D(curves);
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this geometry.
        mini = [inf inf];
        maxi = [-inf -inf];
        for i = 1:length(obj.Curves)
            boxi = boundingBox(obj.Curves{i});
            mini = min(mini, [boxi.XMin boxi.YMin]);
            maxi = max(maxi, [boxi.XMax boxi.YMax]);
        end
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});

        hh = zeros(1, length(obj.Curves));
        for i = 1:length(obj.Curves)
            hh(i) = draw(ax, obj.Curves{i}, varargin{:});
            if ~isempty(style)
                apply(style, hh(i));
            end
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        curves2 = cell(1, length(obj.Curves));
        for i = 1:length(obj.Curves)
            curves2{i} = scale(obj.Curves{i}, factor);
        end
        res = MultiCurve2D(curves2);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        curves2 = cell(1, length(obj.Curves));
        for i = 1:length(obj.Curves)
            curves2{i} = translate(obj.Curves{i}, shift);
        end
        res = MultiCurve2D(curves2{:});
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this geometry.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.

        curves2 = cell(1, length(obj.Curves));
        for i = 1:length(obj.Curves)
            curves2{i} = rotate(obj.Curves{i}, angle, varargin{:});
        end
        res = MultiCurve2D(curves2);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        strArray = cell(1, length(obj.Curves));
        for i = 1:length(obj.Curves)
            strArray{i} = toStruct(obj.Curves{i});
        end
        str = struct('Type', 'MultiCurve2D', 'Curves', {strArray});
    end
end
methods (Static)
    function obj = fromStruct(str)
        % Create a new instance from a structure.
        strArray = str.Curves;
        curveArray = cell(1, length(strArray));
        for i = 1:length(strArray)
            curveArray{i} = Geometry.fromStruct(strArray{i});
        end
        obj = MultiCurve2D(curveArray{:});
    end
end

end % end classdef


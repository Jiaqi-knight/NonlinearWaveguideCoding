classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Point2D < Geometry2D
% A point in the 2-dimensional plane.
%
%   Usage:
%   P = Point2D(X, Y)
%   P = Point2D(COORDS)
%   P = Point2D(PT)
%
%   Example
%   Point2D
%
%   See also
%     Geometry2D, Point3D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    % The X-coordinate of this point.
    X = 0;
    % The Y-coordinate of this point.
    Y = 0;
end 


%% Constructor
methods
    function obj = Point2D(varargin)
        % Constructor for Point2D class.
        
        % empty constructor -> initialize to origin
        if isempty(varargin)
            return;
        end
        
        % copy constructor
        if isa(varargin{1}, 'Point2D')
            that = varargin{1};
            obj.X = that.X;
            obj.Y = that.Y;
            return;
        end
        
        % initialisation constructor with one argument
        if nargin == 1
            var1 = varargin{1};
            if isnumeric(var1) && ~any(size(var1) ~= [1 2])
                obj.X = var1(1);
                obj.Y = var1(2);
            else
                error('Can not parse input for Point2D');
            end
        end
        
        % initialisation from two scalar double argument
        if nargin == 2
            var1 = varargin{1};
            var2 = varargin{2};
            if isnumeric(var1) && isnumeric(var2) && isscalar(var1) && isscalar(var2)
                obj.X = varargin{1};
                obj.Y = varargin{2};
            else
                error('Can not parse inputs for Point2D');
            end
        end
    end

end % end constructors

%% Methods specific to Point2D
methods
    function d = distance(obj, that)
        % Distance between two points.

        % check input
        if ~isa(obj, 'Point2D') || ~isa(that, 'Point2D')
            error('Both arguments must be instances of Point2D');
        end
        
        % compute distance
        dx = obj.X - that.X;
        dy = obj.Y - that.Y;
        d = hypot(dx, dy);
    end
end

%% Methods implementing the Geometry2D interface
methods
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        box = Box2D([obj.X obj.X obj.Y obj.Y]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this point, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'bo'};
        end
        
        % plot line segment
        hh = plot(ax, obj.X, obj.Y, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = Point2D(transformPoint(transform, [obj.X obj.Y]));
    end
end

%% Methods implementing the Geometry2D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = Point2D([obj.X obj.Y] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = Point2D(bsxfun(@plus, [obj.X obj.Y], shift));
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this point.
        %
        % PT2 = rotate(PT, THETA)
        % PT2 = rotate(PT, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.
        
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        % TODO: replace by AffineTransform2D instance
        rot = createRotation(origin, deg2rad(angle));
        verts = transformPoint([obj.X obj.Y], rot);
        
        res = Point2D(verts);
    end
    
    function res = uminus(obj)
        res = Point2D(-obj.X, -obj.Y);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('type', 'Point2D', 'X', obj.X, 'Y', obj.Y);
    end
end
methods (Static)
    function point = fromStruct(str)
        % Create a new instance from a structure.
        point = Point2D(str.X, str.Y);
    end
end


%% sub-indexing methods
methods
    function varargout = subsref(obj, subs)
        % Overrides subsref function for Point2D objects.
        
        % extract reference type
        s1 = subs(1);
        type = s1.type;
        
        % switch between reference types
        if strcmp(type, '.')
            % in case of dot reference, use builtin
            
            % check if we need to return output or not
            if nargout > 0
                % if some output arguments are asked, pre-allocate result
                varargout = cell(nargout, 1);
                [varargout{:}] = builtin('subsref', obj, subs);
            else
                % call parent function, and eventually return answer
                builtin('subsref', obj, subs);
                if exist('ans', 'var')
                    varargout{1} = ans; %#ok<NOANS>
                end
            end
            
        elseif strcmp(type, '()')
            % In case of parens reference, index the inner data
            varargout{1} = 0;
            
            % different processing if 1 or 2 indices are used
            ns = length(s1.subs);
            if ns == 1
                % use index as dimension
                sub1 = s1.subs{1};
                res = zeros(size(sub1));
                res(sub1 == 1) = obj.X;
                res(sub1 == 2) = obj.Y;
                varargout{1} = res;
                
            else
                error('Point2D:subsref', ...
                    'Requires single indices.');
            end
            
        elseif strcmp(type, '{}')
            error('Point2D:subsref', ...
                'can not manage braces reference');
        else
            error('Point2D:subsref', ...
                ['can not manage such reference: ' type]);
        end
    end
end

end % end classdef


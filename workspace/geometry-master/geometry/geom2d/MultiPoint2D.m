classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MultiPoint2D < Geometry2D
% A set of points in the plane.
%
%   Represents a set of Coordinates. 
%
%   Data are stored as a NV-by-2 array.
%
%   Example
%   MultiPoint2D([0 0; 10 0; 10 10]; 0 10]);
%
%   See also
%     Geometry2D, Polygon2D, Point2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRAE - Cepia Software Platform.


%% Properties
properties
    % the point coordinates, given as a N-by-2 array of double
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = MultiPoint2D(varargin)
    % Constructor for MultiPoint2D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 2) ~= 2
                error('Creating a MultiPoint requires an array with two columns, not %d', size(var1, 2));
            end
            obj.Coords = var1;

        else
            obj.Coords = [];
        end
    end

end % end constructors

%% Methods specific to MultiPoint2D
methods
    function centro = centroid(obj)
        % Compute centroid of the points within obj multi-point.
        centro = Point2D(mean(obj.Coords, 1));
    end
end

%% Methods
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = MultiPoint2D(transformPoint(transform, obj.Coords));
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'bo'};
        end
        
        % plot line segment
        xdata = obj.Coords(:,1);
        ydata = obj.Coords(:,2);
        hh = plot(ax, xdata, ydata, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        res = MultiPoint2D(obj.Coords * factor);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        res = MultiPoint2D(bsxfun(@plus, obj.Coords, shift));
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this geometry.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.
        
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        rot = createRotation(origin, deg2rad(angle));
        verts = transformPoint(obj.Coords, rot);
        
        res = MultiPoint2D(verts);
    end
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'MultiPoint2D', 'Coordinates', obj.Coords);
    end
end
methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        if isfield(str, 'Coordinates')
            poly = MultiPoint2D(str.Coordinates);
        elseif isfield(str, 'coordinates')
            poly = MultiPoint2D(str.coordinates);
        else
            error('Field <Coordinates> of MultiPoint2D is not defined');
        end
    end
end


%% sub-indexing methods
methods
    function varargout = size(obj, varargin)
        % Return the size of this multi-point.
        if nargout <= 1
            % compute dim
            if nargin == 2
                s = size(obj.Coords, varargin{1});
            else
                s = size(obj.Coords);
            end
            varargout = {s};
            
        else
            
            s = [size(obj.Coords, 1), size(obj.Coords, 2)];
            varargout = num2cell(s);
        end
    end
    
    function ind = end(obj, k, n)
        % Determine last index when accessing a multi-point.
        %
        %   See Also
        %     subsref
        
        if n == 1
            ind = size(obj.Coords, 1);
        elseif n == 2
            if k == 1 
                ind = size(obj.Coords, 1);
            elseif k == 2
                ind = 2;
            else
                error('MultiPoint2D:end', ...
                    'Only two dimensions can be accessed');
            end
        else
            error('MultiPoint2D:end', ...
                'Only two dimensions can be accessed');
        end
    end
    
    function varargout = subsref(obj, subs)
        % Overrides subsref function for LineString2D objects.
        
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
                % Return the point(s) at specified index(ices)
                sub1 = s1.subs{1};
                if length(sub1) == 1
                    varargout{1} = Point2D(obj.Coords(sub1, :));
                else
                    varargout{1} = MultiPoint2D(obj.Coords(sub1, :));
                end
                
            elseif ns == 2
                % Return coordinates at specified indices and dimensions
                sub1 = s1.subs{1};
                if ischar(sub1) && strcmp(sub1, ':')
                    sub1 = 1:size(obj.Coords, 1);
                end
                sub2 = s1.subs{2};
                if ischar(sub2) && strcmp(sub2, ':')
                    sub2 = 1:2;
                end
                varargout{1} = obj.Coords(sub1, sub2);
                
            else
                error('MultiPoint2D:subsref', ...
                    'Too many indices.');
            end
            
        elseif strcmp(type, '{}')
            error('MultiPoint2D:subsref', ...
                'can not manage braces reference');
        else
            error('MultiPoint2D:subsref', ...
                ['can not manage such reference: ' type]);
        end
    end
end

end % end classdef


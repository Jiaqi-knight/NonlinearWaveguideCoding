classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MultiPoint3D < handle
% A set of points in the 3D space.
%
%   Class MultiPoint3D
%
%   Example
%   MultiPoint3D
%
%   See also
%     Point3D, MultiPoint2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-07,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % the vertex coordinates, given as a N-by-3 array of double
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = MultiPoint3D(varargin)
    % Constructor for MultiPoint3D class.

        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 2) ~= 3
                error('Creating a 3D MultiPoint requires an array with three columns, not %d', size(var1, 2));
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
        centro = Point3D(mean(obj.Coords, 1));
    end
end


%% Methods
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = MutliPoint3D(transformPoint(transform, obj.Coords));
    end
    
    function box = boundingBox(obj)
        % Returns the bounding box of this shape.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % draw the geometric primitive
        hh = plot3(ax, obj.Coords(:,1), obj.Coords(:,2), obj.Coords(:,3), varargin{:});

        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
        
        % format output argument
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = MultiPoint3D(obj.Coords * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = MultiPoint3D(bsxfun(@plus, obj.Coords, shift));
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'MultiPoint3D', 'Coordinates', obj.Coords);
    end
end
methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        poly = MultiPoint3D(str.Coordinates);
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
                ind = 3;
            else
                error('MultiPoint3D:end', ...
                    'Only two dimensions can be accessed');
            end
        else
            error('MultiPoint3D:end', ...
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
                    varargout{1} = Point3D(obj.Coords(sub1, :));
                else
                    varargout{1} = MultiPoint3D(obj.Coords(sub1, :));
                end
                
            elseif ns == 2
                % Return coordinates at specified indices and dimensions
                sub1 = s1.subs{1};
                if ischar(sub1) && strcmp(sub1, ':')
                    sub1 = 1:size(obj.Coords, 1);
                end
                sub2 = s1.subs{2};
                if ischar(sub2) && strcmp(sub2, ':')
                    sub2 = 1:3;
                end
                varargout{1} = obj.Coords(sub1, sub2);
                
            else
                error('MultiPoint3D:subsref', ...
                    'Too many indices.');
            end
            
        elseif strcmp(type, '{}')
            error('MultiPoint3D:subsref', ...
                'can not manage braces reference');
        else
            error('MultiPoint3D:subsref', ...
                ['can not manage such reference: ' type]);
        end
        
    end
end

end % end classdef


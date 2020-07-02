classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LinearRing3D < Geometry3D
% A 3D closed polyline composed of several line segments.
%
%   Represents a polyline defined be a series of Coords. 
%
%   Data are represented by a NV-by-3 array.
%
%   Example
%      ring = LinearRing3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);
%      draw(ring);
%
%   See also
%     Geometry3D, LineString3D, LinearRing2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-09-11,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Properties
properties
    % the set of Coords, given as a N-by-3 array of coordinates
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = LinearRing3D(varargin)
    % Constructor for LineString3D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 2) ~= 3
                error('Creating a LinearRing3D requires an array with three columns, not %d', size(var1, 2));
            end
            obj.Coords = var1;

        else
            obj.Coords = [];
        end
    end

end % end constructors


%% Methods specific to LinearRing3D
methods
    function res = smooth(obj, M)
        %SMOOTH Smooth a linear ring using local averaging.
        %
        %   RING2 = smooth(RING, M)

        % create convolution vector
        v2 = ones(M, 1) / M;
        
        % compute the number of elements before and after
        M1 = floor((M - 1) / 2);
        M2 = ceil((M - 1) / 2);
        
        % allocate memory for result
        res = zeros(size(obj.Coords));
        
        % iterate over dimensions
        for d = 1:3
            vals = [obj.Coords(end-M1+1:end,d); obj.Coords(:,d) ; obj.Coords(1:M2,d)];
            resd = conv(vals, v2, 'same');
            res(:,d) = resd(M1+1:end-M2);
        end

        % convert result to LineString object
        res = LinearRing3D(res);
    end
    
    function l = length(obj)
        % Compute the length of this polyline.
        %
        % L = length(obj);

        % compute the sum of the length of each line segment
        l = sum(sqrt(sum(diff(obj.Coords([1:end 1], :)).^2, 2)));
    end
    
    function centro = centroid(obj)
        % Compute the curve centroid of this polyline.
        
        % compute center and length of each line segment
        centers = (obj.Coords(1:end,:) + obj.Coords([2:end 1],:))/2;
        lengths = sqrt(sum(diff(obj.Coords([1:end 1], :)).^2, 2));
        
        % centroid of edge centers weighted by edge lengths
        centro = Point3D(sum(bsxfun(@times, centers, lengths), 1) / sum(lengths));
    end
end

%% Data access specific to Linear Ring 
methods
    function verts = vertices(obj)
        % Return vertices as a new instance of MultiPoint3D.
        verts = MultiPoint3D(obj.Coords);
    end
    
    function nv = vertexNumber(obj)
        % Get the number of vertices in this polyline.
        nv = size(obj.Coords, 1);
    end
end



%% Methods generic to curve objects
methods
    function res = reverse(obj)
        res = LinearRing3D(obj.Coords([1 end:-1:2],:));
    end
end


%% Methods
methods
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function varargout = draw(varargin)
        % DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        inds = [1:size(obj.Coords, 1) 1]; 
        h = plot3(ax, obj.Coords(inds,1), obj.Coords(inds,2), obj.Coords(inds,3), varargin{:});
        
        if ~isempty(style)
            apply(style, h);
        end
        
        if nargout > 0
            varargout = {h};
        end
     end
    
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = LinearRing3D(transformPoint(transform, obj.Coords));
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = LinearRing3D(obj.Coords * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = LinearRing3D(bsxfun(@plus, obj.Coords, shift));
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LinearRing3D', 'Coordinates', obj.Coords);
    end
end

methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        poly = LinearRing3D(str.Coordinates);
    end
end


%% sub-indexing methods
methods
    function varargout = size(obj, varargin)
        % Return the size of this polyline.
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
        % Determine last index when accessing a polyline.
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
                error('LinearRing3D:end', ...
                    'Only two dimensions can be accessed');
            end
        else
            error('LinearRing3D:end', ...
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
                error('LinearRing3D:subsref', ...
                    'Too many indices.');
            end
            
        elseif strcmp(type, '{}')
            error('LinearRing3D:subsref', ...
                'can not manage braces reference');
        else
            error('LinearRing3D:subsref', ...
                ['can not manage such reference: ' type]);
        end
    end
end

end % end classdef


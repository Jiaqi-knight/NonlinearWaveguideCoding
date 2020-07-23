classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LinearRing2D < Curve2D
% A closed polyline in the plane.
%
%   Represents a linear ring defined be a series of Coords. 
%
%   Data are represented by a NV-by-2 array.
%
%   Example
%     ring = LinearRing2D([0 0; 10 0; 10 10; 0 10]);
%     draw(ring);
%
%   See also
%     Geometry2D, LineString2D, Polygon2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Properties
properties
    % the set of vertex coordinates, given as a N-by-2 array of double
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = LinearRing2D(varargin)
    % Constructor for LinearRing2D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 2) ~= 2
                error('Creating a linear ring requires an array with two columns, not %d', size(var1, 2));
            end
            obj.Coords = var1;

        else
            obj.Coords = [];
        end
    end

end % end constructors

%% Methods specific to LinearRing2D
methods
    function res = simplify(obj, tol)
        % Douglas-Peucker simplification of this linear ring.
        %
        %   LS2 = simplify(LD, TOL);
        %
        %   Example
        %     elli = Ellipse2D([50 50 40 30  20]);
        %     ring = asPolyline(elli, 200);
        %     figure; axis equal; axis([0 100 0 100]); hold on;
        %     draw(ring, 'b');
        %     ring2 = simplify(ring, 2);
        %     draw(ring2, 'm');
        %     drawVertices(ring2);
        %
        %   See also
        %     resample
        
        % Use LineString2D implementation
        ls = LineString2D(obj.Coords([1:end 1], :));
        ls2 = simplify(ls, tol);
        res = LinearRing2D(ls2.Coords(1:end-1, :));
    end
    
    function poly2 = resample(obj, n)
        % RESAMPLE Resample this polyline with a given number of vertices.
        %
        %   Syntax:  POLY2 = resample(POLY, N);

        % compute arc length along each ertex
        s = verticesArcLength(obj);
        
        % distribute N+1 points equally spaced (the last one is removed at
        % the end)
        Lmax = s(end);
        pos = linspace(0, Lmax, n+1);

        coords = obj.Coords([1:end 1], :);
        poly2 = zeros(n+1, size(coords, 2));
        for i = 1:n+1
            % index of surrounding vertices before and after
            ind0 = find(s <= pos(i), 1, 'last');
            ind1 = find(s >= pos(i), 1, 'first');
            
            if ind0 == ind1
                % get position of a vertex in input polyline
                poly2(i, :) = coords(ind0, :);
                continue;
            end
            
            % position of surrounding vertices
            pt0 = coords(ind0, :);
            pt1 = coords(ind1, :);
            
            % weights associated to each neighbor
            l0 = pos(i) - s(ind0);
            l1 = s(ind1) - pos(i);
            
            % linear interpolation of neighbor positions
            if (l0 + l1) > Lmax * 1e-12
                poly2(i, :) = (pt0 * l1 + pt1 * l0) / (l0 + l1);
            else
                % if neighbors are too close, do not use interpolation
                poly2(i, :) = pt0;
            end
        end
        
        % Remove last vertex (same as first one) and convert result to a
        % LinearRing2D instance 
        poly2 = LinearRing2D(poly2(1:end-1,:));
    end
    
    function al = verticesArcLength(obj)
        % Return the arc length at each vertex of the polyline.
        %
        % Syntax:
        %   AL = verticesArcLength(POLY);
        % Returns an array AL with one value more than the number of
        % vertices, as the first vertex is counted twice.
        
        % compute the cumulative  sum of the length of each line segment,
        % and add 0 for the first vertex.
        al = [0 ; cumsum(sqrt(sum(diff(obj.Coords([1:end 1], :)).^2, 2)))];
    end
    
    function centro = centroid(obj)
        % Return the centroid of this polyline.
        
        % compute center and length of each line segment
        centers = (obj.Coords(1:end,:) + obj.Coords([2:end 1],:))/2;
        lengths = sqrt(sum(diff(obj.Coords([1:end 1], :), 1).^2, 2));
        
        % centroid of edge centers weighted by edge lengths
        centro = Point2D(sum(bsxfun(@times, centers, lengths), 1) / sum(lengths));
    end
end


%% Methods or management of vertices
methods
    function verts = vertices(obj)
        % Return vertices as a new instance of MultiPoint2D.
        verts = MultiPoint2D(obj.Coords);
    end
    
    function nv = vertexNumber(obj)
        % Get the number of vertices.
        nv = size(obj.Coords, 1);
    end
end


%% Methods generic to curve objects
methods
    function p = length(obj)
        % Compute the length of this polyline.
        dp = diff(obj.Coords([1:end 1], :), 1, 1);
        p = sum(hypot(dp(:, 1), dp(:, 2)));
    end
    
    function res = reverse(obj)
        res = LinearRing2D(obj.Coords([1 end:-1:2],:));
    end
end


%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = LinearRing2D(transformPoint(transform, obj.Coords));
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function h = drawVertices(varargin)
        % Draw vertices of this polyline, with optional drawing options.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');
        
        % default options
        if isempty(varargin)
            varargin = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'w'};
        end
        
        % extract data
        xdata = obj.Coords(:,1);
        ydata = obj.Coords(:,2);
        
        hh = plot(ax, xdata, ydata, varargin{:});
        if ~isempty(style)
            apply(style, hh);
        end
        
        if holdState
            hold(ax, 'on');
        else
            hold(ax, 'off');
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');

        % default options
        drawLines = true;
        drawVertices = false;
        if ~isempty(style)
            drawLines = style.LineVisible;
            drawVertices = style.MarkerVisible;
        end
        
        % parse some options
        inds = strcmpi(varargin, 'drawVertices');
        if any(inds)
            inds = find(inds(1));
            drawVertices = varargin{inds+1};
            varargin([inds inds+1]) = [];
        end
        
        % extract data
        xdata = obj.Coords(:,1);
        ydata = obj.Coords(:,2);
        
        % draw outline
        h1 = [];
        if drawLines
            if isempty(varargin)
                varargin = {'Color', 'b', 'LineStyle', '-'};
            end
            h1 = plot(ax, xdata([1:end 1]), ydata([1:end 1]), varargin{:});
            if ~isempty(style)
                apply(style, h1);
            end
        end
        
        % optionnally draw markers
        h2 = [];
        if drawVertices
            options = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none'};
            h2 = plot(ax, xdata, ydata, options{:});
            if ~isempty(style)
                apply(style, h2);
            end
        end
        
        if holdState
            hold(ax, 'on');
        else
            hold(ax, 'off');
        end
        
        if nargout > 0
            h = [h1 h2];
        end
    end
end

%% Methods implementing the Geometry2D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = LinearRing2D(obj.Coords * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = LinearRing2D(bsxfun(@plus, obj.Coords, shift));
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this polyline.
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
        
        res = LinearRing2D(verts);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LinearRing2D', 'Coordinates', obj.Coords);
    end
end
methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        if isfield(str, 'Coordinates')
            poly = LinearRing2D(str.Coordinates);
        elseif isfield(str, 'coordinates')
            poly = LinearRing2D(str.coordinates);
        else
            error('Field <Coordinates> of LinearRing2D is not defined');
        end
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
                ind = 2;
            else
                error('LinearRing2D:end', ...
                    'Only two dimensions can be accessed');
            end
        else
            error('LinearRing2D:end', ...
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
                error('LinearRing2D:subsref', ...
                    'Too many indices.');
            end
            
        elseif strcmp(type, '{}')
            error('LinearRing2D:subsref', ...
                'can not manage braces reference');
        else
            error('LinearRing2D:subsref', ...
                ['can not manage such reference: ' type]);
        end
    end
end

end % end classdef


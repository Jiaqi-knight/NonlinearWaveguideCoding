classdef StraightLine2D < Geometry2D
% A straight line, unbounded in each direction.
%
%   Class StraightLine2D
%
%   Example
%     P1 = Point2D([40 10]);
%     P1 = Point2D([20 40]);
%     L = StraightLine2D(P1, P2);
%     figure; axis equal;axis([0 50 0 50]; holdon;
%     draw(L);
%
%   See also
%     LineSegment2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-22,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
    % The origin of this straight line, as a 1-by-2 numeric array.
    Origin = [0 0];
    
    % The direction vector of this straight line, as a 1-by-2 numeric array.
    Direction = [1 0];
    
end % end properties


%% Constructor
methods
    function obj = StraightLine2D(varargin)
        % Constructor for StraightLine2D class
        
        if nargin == 0
            % empty constructor -> already initialized to origin
            
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'StraightLine2D')
                % copy constructor
                var1 = varargin{1};
                obj.Origin    = var1.Origin;
                obj.Direction = var1.Direction;
            elseif isnumeric(var1) && all(size(var1) == [1 4])
                % numeric input consider [X0 Y0 DX DY].
                obj.Origin    = [var1(1) var1(2)];
                obj.Direction = [var1(3) var1(4)];
            else
                error('Can not parse input for StraightLine2D');
            end
            
        elseif nargin == 2
            % initialisation from two arguments
            
            % parse origin
            var1 = varargin{1};
            if isa(var1, 'Point2D')
                obj.Origin = [var1.X var1.Y];
            elseif isnumeric(var1)
                obj.Origin = var1(1, 1:2);
            else
                error('Can not interpret second argument');
            end
            
            % second argument can be either another point, or the direction
            var2 = varargin{2};
            if isa(var2, 'Point2D')
                obj.Direction = [var2.X var2.Y] - obj.Origin;
            elseif isa(var2, 'Vector3D')
                obj.Direction = [var2.X var2.Y];
            elseif isnumeric(var2)
                % numeric inpu consider another point as default.
                obj.Direction = var2 - obj.Origin;
            else
                error('Can not interpret second argument');
            end
            
        else
            error('Wrong number of input arguments.');
        end
        
    end

end % end constructors

%% Methods
methods
    function p = origin(obj)
        p = Point2D(obj.Origin);
    end
    
    function v = direction(obj)
        v = Vector2D(obj.Direction);
    end
end

%% Methods
methods
    function edge = clip(obj, box)
        % Compute the portion of line visible within a box.
        
        % extract limits of the box
        xmin = box.XMin;
        xmax = box.XMax;
        ymin = box.YMin;
        ymax = box.YMax;
        
        % use direction vector for box edges similar to direction vector of the
        % line in order to reduce computation errors
        delta = norm(direction(obj));
        
        % compute intersection with each edge of the box
        px1 = intersection(obj, StraightLine2D([xmin ymin  delta 0]));  % lower edge
        px2 = intersection(obj, StraightLine2D([xmax ymin 0  delta]));  % right edge
        py1 = intersection(obj, StraightLine2D([xmax ymax -delta 0]));  % upper edge
        py2 = intersection(obj, StraightLine2D([xmin ymax 0 -delta]));  % left edge
        
        % remove undefined intersections (case of lines parallel to box edges)
        points = {};
        if ~isempty(px1), points = [points {px1}]; end
        if ~isempty(px2), points = [points {px2}]; end
        if ~isempty(py1), points = [points {py1}]; end
        if ~isempty(py2), points = [points {py2}]; end
        
        % sort points according to their position on the line
        pos = zeros(1, length(points));
        for i = 1:length(points)
            pos(i) = position(obj, points{i});
        end
        [pos, inds] = sort(pos); %#ok<ASGLU>
        points = points(inds);
        
        % create clipped edge by using the two points in the middle
        ind = length(points) / 2;
        inter1 = points{ind};
        inter2 = points{ind+1};
        edge = LineSegment2D(inter1, inter2);
        
        % check that middle point of the edge is contained in the box
        mid = middlePoint(edge);
        xOk = xmin <= mid.X && mid.X <= xmax;
        yOk = ymin <= mid.Y && mid.Y <= ymax;
        
        % if one of the bounding condition is not met, set edge to empty
        if ~(xOk && yOk)
            edge = [];
        end
    end
end % end methods

methods
    function p = intersection(line1, line2)
        % Compute the intersection point of two lines.
        
        tol = 1e-12;
        
        % Check lines are not parallel
        denom = line1.Direction(1) .* line2.Direction(2) - line2.Direction(1) .* line1.Direction(2);
        if abs(denom) < tol
            p = [];
            return;
        end
        
        % Extract coordinates of intersecting lines
        x1 =  line1.Origin(1);
        y1 =  line1.Origin(2);
        dx1 = line1.Direction(1);
        dy1 = line1.Direction(2);
        
        % extract base coordinates of second lines
        x2 =  line2.Origin(1);
        y2 =  line2.Origin(2);
        dx2 = line2.Direction(1);
        dy2 = line2.Direction(2);
        
        % compute coordinate differences of origin points
        dx = x2 - x1;
        dy = y2 - y1;
        
        % Compute coordinates of intersection point
        x0 = (x2 * dy2 * dx1 - dy * dx1 * dx2 - x1 * dy1 * dx2) / denom ;
        y0 = (dx * dy1 * dy2 + y1 * dx1 * dy2 - y2 * dx2 * dy1) / denom ;
        
        % concatenate result
        p = Point2D(x0, y0);
    end
    
    function pos = position(obj, point)
        % Compute position of a point on this line.
        
        % direction vector of the lines
        vx = obj.Direction(1);
        vy = obj.Direction(2);
        
        % difference of coordinates between point and line origins
        dx = point.X - obj.Origin(1);
        dy = point.Y - obj.Origin(2);
        
        % squared norm of direction vector, with a check of validity
        delta = vx .* vx + vy .* vy;
        if delta < eps
            pos = 0;
            return
        end
        
        % compute position of points projected on the line, by using normalised dot
        % product (NP-by-NL array)
        pos = (dx * vx + dy * vy) / delta;
    end
end


%% Methods implementing the Geometry2D interface
methods
    function box = boundingBox(obj) %#ok<STOUT,MANU>
        % Return the bounding box of this geometry.
        error('Straight line is not a bounded geometry.');
    end
    
    function h = draw(varargin)
        %DRAW Draw this point, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default style for drawing lines
        if length(varargin) ~= 1
            varargin = [{'color', 'b'}, varargin];
        end

        % extract bounding box of the current axis
        xlim = get(ax, 'xlim');
        ylim = get(ax, 'ylim');
        
        % clip lines with current axis box
        edge = clip(obj, [xlim ylim]);
        
        if ~isempty(ege)
            % draw current clipped line
            hh = draw(ax, edge, varargin{:});
        
            if ~isempty(style)
                apply(style, hh);
            end
        else
            hh = -1;
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        origin = transformPoint(transform, obj.Origin);
        direction = transformVector(transform, obj.Direction);
        res = StraightLine2D(origin, direction);
    end
end


%% Methods implementing the Geometry2D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        if ~isscalar(factor)
            error('Requires scaling factor to be a scalar');
        end
        origin = obj.Origin * factor;
        direction = obj.Direction * factor;
        res = StraightLine2D(origin, direction);
    end
    
    function res = rotate(obj, varargin)
        % Return a rotated version of this geometry.
        p1 = rotate(origin(obj), varargin{:});
        p2 = rotate(Point2D(obj.Origin + obj.Direction), varargin{:});
        res = StraightLine2D(p1, p2);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        origin = obj.Origin + shift;
        res = StraightLine2D(origin, obj.Direction);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'StraightLine2D', ...
            'Origin', obj.Origin, ...
            'Direction', obj.Direction);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        line = StraightLine2D(str.Origin, str.Origin+str.Direction);
    end
end

end % end classdef

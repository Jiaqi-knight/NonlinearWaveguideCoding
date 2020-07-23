classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) TransformedGrid2D < Geometry2D
%TRANSFORMEDGRID2D Utility class for representing deformation of 2D grids.
%
%   Class TransformedGrid2D
%
%   Example
%     % Display a square grid rotated around center of image
%     % (uses AffineTransform class from MatRegister library)
%     img = imread('cameraman.tif');
%     figure; imshow(img); hold on;
%     tra = AffineTransform.createTranslation([128 128]);
%     rot = AffineTransform.createRotation(pi/6);
%     transfo = tra * rot * inverse(tra);
%     tg = TransformedGrid2D(0:32:256, 0:32:256, transfo, [8 8]);
%     draw(tg, 'b')
%
%   See also
%     TransformedGrid3D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The x-coordinates of the vertices, as a 1-by-M array.
    XVertices;
    % The y-coordinates of the vertices, as a 1-by-N array.
    YVertices;
    
    % The transform to apply to the grid.
    Transform;
    
    % The discretization step of grid lines in the X direction.
    DX;
    % The discretization step of grid lines in the Y direction.
    DY;
    
end % end properties


%% Constructor
methods
    function obj = TransformedGrid2D(LX, LY, transfo, steps)
        % Constructor for TransformedGrid2D class.
        %
        % OBJ = TransformedGrid2D(LX, LY, TRANSFO, STEPS)
        % OBJ = TransformedGrid2D(LX, LY, TRANSFO)
        % Automatically determines steps.
        %
        
        obj.XVertices = LX;
        obj.YVertices = LY;
        obj.Transform = transfo;
        if nargin > 3
            obj.DX = steps(1);
            obj.DY = steps(2);
        else
            obj.DX = (LX(2)-LX(1)) / 10;
            obj.DY = (LY(2)-LY(1)) / 10;
        end            
    
    end

end % end constructors


%% Methods implementing the Geometry interface
methods
    function res = transform(obj, transform) %#ok<STOUT>
        % Apply a geometric transform to this geometry.
        error('Transform not implemented for TransformedGrid objects');
    end
    
    function box = boundingBox(obj)
        % Returns the bounding box of this geometry.
        [x, y] = meshgrid(obj.XVertices, obj.YVertices);
        coordsT = transformPoint(obj.Transform, [x(:) y(:)]);
        x2 = coordsT(:,1);
        y2 = coordsT(:,2);
        box = Box2D([min(x2(:)) max(x2(:)) min(y2(:)) max(y2(:))]);
    end

    function h = draw(varargin)
        % Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'color', [.5 .5 .5], 'linestyle', '-'};
        end
        
        % grille horizontale
        lx = obj.XVertices(1):obj.DX:obj.XVertices(end);
        ly = obj.YVertices;
        hy = zeros(1,length(ly));
        for iy = 1:length(ly)
            coords = [lx(:) ones(size(lx(:)))*ly(iy)];
            coords = transformPoint(obj.Transform, coords);
            poly = LineString2D(coords);
            hy(iy) = draw(ax, poly, varargin{:});
            if ~isempty(style)
                apply(style, hy(iy));
            end
        end
        
        % grille verticale
        lx = obj.XVertices;
        ly = obj.YVertices(1):obj.DY:obj.YVertices(end);
        hx = zeros(1,length(lx));
        for ix = 1:length(lx)
            coords = [ones(size(ly(:)))*lx(ix) ly(:)];
            coords = transformPoint(obj.Transform, coords);
            poly = LineString2D(coords);
            hx(ix) = draw(ax, poly, varargin{:});
            if ~isempty(style)
                apply(style, hx(ix));
            end
        end
        
        if nargout > 0
            h = [hx hy];
        end
    end
 
    
    function res = scale(obj, factor)
        % Returns a scaled version of this geometry.
        vx2 = obj.XVertices * factor;
        vy2 = obj.YVertices * factor;
        steps2 = [obj.DX * factor, obj.DY * factor];
        res = TransformedGrid2D(vx2, vy2, obj.Transform, steps2);
    end
    
    function res = translate(obj, shift)
        % Returns a translated version of this geometry.
        vx2 = obj.XVertices + shift(1);
        vy2 = obj.YVertices + shift(2);
        steps2 = [obj.DX obj.DY];
        res = TransformedGrid2D(vx2, vy2, obj.Transform, steps2);
    end
    
    function rotate(obj, angle, varargin)
        % Non supported operation.
        error('Non supported operation');
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'TransformedGrid2D', ...
            'XVertices', obj.XVertices, ...
            'YVertices', obj.YVertices, ...
            'Transform', obj.Transform, ...
            'DX', obj.DX, ...
            'DY', obj.DY);
    end
end
methods (Static)
    function circ = fromStruct(str)
        % Create a new instance from a structure.
        transfo = Transform.fromStruct(str);
        circ = TransformedGrid2D(str.XVertices, str.YVertices, transfo, [str.DX str.DY]);
    end
end


end % end classdef


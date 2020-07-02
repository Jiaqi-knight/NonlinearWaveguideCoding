classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) TransformedGrid3D < Geometry3D
%TRANSFORMEDGRID3D Utility class for representing deformation of 3D grids.
%
%   Class TransformedGrid3D
%
%   Example
%     TransformedGrid3D
%
%   See also
%     Transform, TransformedGrid2D

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
    % The z-coordinates of the vertices, as a 1-by-P array.
    ZVertices;
    
    % The transform to apply to the grid.
    Transform;
    
    % The discretization step of grid lines in the X direction.
    DX;
    % The discretization step of grid lines in the Y direction.
    DY;
    % The discretization step of grid lines in the Z direction.
    DZ;
    
end % end properties


%% Constructor
methods
    function obj = TransformedGrid3D(LX, LY, LZ, transfo, steps)
        % Constructor for TransformedGrid3D class.
        %
        % OBJ = TransformedGrid2D(LX, LY, LZ, TRANSFO, STEPS)
        % Create a new instance for a 3D M-by-N-by-P grid transformed by
        % the transform TRANSFO.
        % LX, LY and LZ are row vectors specifying grid coordinates in
        % original basis. TRANSFO is an instance of Transform. STEPS is a
        % 1-by-3 row vector specifying discretization step of grid lines.
        %
        % OBJ = TransformedGrid2D(LX, LY, LZ, TRANSFO)
        % Automatically determines steps.
        %
        
        obj.XVertices = LX;
        obj.YVertices = LY;
        obj.ZVertices = LZ;
        obj.Transform = transfo;
        
        if nargin > 4
            obj.DX = steps(1);
            obj.DY = steps(2);
            obj.DZ = steps(3);
        else
            obj.DX = (LX(2)-LX(1)) / 10;
            obj.DY = (LY(2)-LY(1)) / 10;
            obj.DZ = (LZ(2)-LZ(1)) / 10;
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
        % Return the bounding box of this geometry.
        [x, y, z] = meshgrid(obj.XVertices, obj.YVertices, obj.ZVertices);
        box = Box3D([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
    end

    function h = draw(varargin)
        % Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'color', [.5 .5 .5], 'linestyle', '-'};
        end
        
        % grid in the X direction
        lx = obj.XVertices(1):obj.DX:obj.XVertices(end);
        ly = obj.YVertices;
        lz = obj.ZVertices;
        
        tmp = ones(size(lx(:)));
        hx = zeros(1, length(ly) * length(lz));
        
        ix = 1;
        for iz = 1:length(lz)
            for iy = 1:length(ly)
                coords = [lx(:), tmp * ly(iy), tmp * lz(iz)];
                coords = transformPoint(obj.Transform, coords);
                
                poly = LineString3D(coords);
                
                hx(ix) = draw(ax, poly, varargin{:});
                if ~isempty(style)
                    apply(style, hx(ix));
                end
                
                ix = ix + 1;
            end
        end
        
        % grid in the Y direction
        lx = obj.XVertices;
        ly = obj.YVertices(1):obj.DY:obj.YVertices(end);
        lz = obj.ZVertices;
        
        tmp = ones(size(ly(:)));
        hy = zeros(1, length(lx) * length(lz));
        
        iy = 1;
        for iz = 1:length(lz)
            for ix = 1:length(lx)
                coords = [tmp * lx(ix), ly(:), tmp * lz(iz)];
                coords = transformPoint(obj.Transform, coords);
                
                poly = LineString3D(coords);
                
                hy(iy) = draw(ax, poly, varargin{:});
                if ~isempty(style)
                    apply(style, hy(iy));
                end
                
                iy = iy + 1;
            end
        end
        
        % grid in the Z direction
        lx = obj.XVertices;
        ly = obj.YVertices;
        lz = obj.ZVertices(1):obj.DZ:obj.ZVertices(end);
        
        tmp = ones(size(lz(:)));
        hz = zeros(1, length(lx) * length(ly));
        
        iz = 1;
        for iy = 1:length(ly)
            for ix = 1:length(lx)
                coords = [tmp * lx(ix), tmp * ly(iy), lz(:)];
                coords = transformPoint(obj.Transform, coords);
                
                poly = LineString3D(coords);
                
                hz(iz) = draw(ax, poly, varargin{:});
                if ~isempty(style)
                    apply(style, hz(iz));
                end
                
                iz = iz + 1;
            end
        end
           
        if nargout > 0
            h = [hx, hx, hz];
        end
    end
 
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        vx2 = obj.XVertices * factor;
        vy2 = obj.YVertices * factor;
        vz2 = obj.ZVertices * factor;
        steps2 = [obj.DX * factor, obj.DY * factor, obj.DZ * factor];
        res = TransformedGrid3D(vx2, vy2, vz2, obj.Transform, steps2);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        vx2 = obj.XVertices + shift(1);
        vy2 = obj.YVertices + shift(2);
        vz2 = obj.ZVertices + shift(3);
        steps2 = [obj.DX obj.DY obj.DZ];
        res = TransformedGrid3D(vx2, vy2, vz2, obj.Transform, steps2);
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'TransformedGrid3D', ...
            'XVertices', obj.XVertices, ...
            'YVertices', obj.YVertices, ...
            'ZVertices', obj.ZVertices, ...
            'Transform', obj.Transform, ...
            'DX', obj.DX, ...
            'DY', obj.DY, ...
            'DZ', obj.DZ);
    end
end
methods (Static)
    function tg = fromStruct(str)
        % Create a new instance from a structure
        transfo = Transform.fromStruct(str);
        tg = TransformedGrid3D(str.XVertices, str.YVertices, str.ZVertices, transfo, [str.DX str.DY str.DZ]);
    end
end


end % end classdef


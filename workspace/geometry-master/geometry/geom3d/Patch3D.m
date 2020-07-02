classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Patch3D < Geometry3D
% A 3D parametric surface defined by three arrays x, y, and z.
%
%   Class Patch3D
%
%   Example
%   Patch3D
%
%   See also
%     Patch2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-25,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    X;
    Y;
    Z;
end % end properties


%% Constructor
methods
    function obj = Patch3D(varargin)
    % Constructor for Patch3D class

        if nargin == 3
            obj.X = varargin{1};
            obj.Y = varargin{2};
            obj.Z = varargin{3};
        end
    end

end % end constructors


%% Methods specific to Patch3D
methods
    function res = smooth(obj, M)
        % Smooth this 3D patch using local averaging.

        % create convolution vector
        v2 = ones(M, 1) / M;
        
        X2 = padarray(obj.X, [M M], 'replicate');
        X2 = conv2(X2, v2, 'same');
        X2 = X2(M+1:end-M, M+1:end-M);
        
        Y2 = padarray(obj.Y, [M M], 'replicate');
        Y2 = conv2(Y2, v2, 'same');
        Y2 = Y2(M+1:end-M, M+1:end-M);

        Z2 = padarray(obj.Z, [M M], 'replicate');
        Z2 = conv2(Z2, v2, 'same');
        Z2 = Z2(M+1:end-M, M+1:end-M);

        % convert result to Patch3D object
        res = Patch3D(X2, Y2, Z2);
    end

    function verts = vertices(obj)
        % Return vertices as a new instance of MultiPoint3D.
        coords = [obj.X(:) obj.Y(:) obj.Z(:)];
        verts = MultiPoint3D(coords);
    end
    
    function drawSubGrid(varargin)
        % Draw a grid within this patch.
        %
        %   drawSubGrid(OBJ, 1) simply displays the boundary of the patch.
        %
        
        % isolates the object instance
        ind = find(cellfun(@(x) isa(x, 'Patch3D'), varargin));
        obj = varargin{ind(1)};
        varargin(ind(1)) = [];
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % determines the number of patch tiles
        nTiles = 1;
        if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
            nTiles = varargin{1};
            varargin(1) = [];
        end
        
        % compute indices of curves
        subI = round(linspace(1, size(obj.X,1), nTiles+1));
        subJ = round(linspace(1, size(obj.X,2), nTiles+1));

        % draw 3D curves in first direction
        for i = 1:length(subI)
            x = obj.X(subI(i), :)';
            y = obj.Y(subI(i), :)';
            z = obj.Z(subI(i), :)';

            plot3(hAx, x, y, z, varargin{:});
        end
        
        % draw 3D curves in second direction
        for i = 1:length(subJ)
            x = obj.X(:, subJ(i));
            y = obj.Y(:, subJ(i));
            z = obj.Z(:, subJ(i));
            plot3(hAx, x, y, z, varargin{:});
        end
    end
end

%% Methods implementing the Geometry3D interface
methods
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        xmin = min(obj.X(:));
        xmax = max(obj.X(:));
        ymin = min(obj.Y(:));
        ymax = max(obj.Y(:));
        zmin = min(obj.Z(:));
        zmax = max(obj.Z(:));
        box = Box3D([xmin xmax ymin ymax zmin zmax]);
    end    

    function h = draw(varargin)
        %DRAW Draw the patch, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        hh = surf(ax, 'XData', obj.X, 'YData', obj.Y, 'ZData', obj.Z, varargin{:});

        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
        
        % format output argument
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transfo)
        % Return a transformed version of this geometry.
        if isa(transfo, 'AffineTransform2D')
            coords = [obj.X(:) obj.Y(:)];
            coords = transformPoint(transfo, coords);
            dims = size(obj.X);
            res = Patch3D(reshape(coords(:,1), dims), reshape(coords(:,2), dims), obj.Z);
        else
            coords = [obj.X(:) obj.Y(:) obj.Z(:)];
            coords = transformPoint(transfo, coords);
            dims = size(obj.X);
            res = Patch3D(reshape(coords(:,1), dims), reshape(coords(:,2), dims), reshape(coords(:,3), dims));
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        res = Patch3D(obj.X * factor, obj.Y * factor, obj.Z * factor);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        res = Patch3D(obj.X + shift(1), obj.Y + shift(2), obj.Z + shift(3));
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Patch3D', 'X', obj.X, 'Y', obj.Y, 'Z', obj.Z);
    end
end
methods (Static)
    function obj = fromStruct(str)
        % Create a new instance from a structure.
        obj = Patch3D(str.X, str.Y, str.Z);
    end
end

end % end classdef


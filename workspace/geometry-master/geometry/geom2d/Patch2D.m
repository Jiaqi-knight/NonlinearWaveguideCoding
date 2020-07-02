classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Patch2D < Geometry2D
% A 2D parametric grid defined by two arrays x and y.
%
%   Class Patch2D
%
%   Example
%   Patch2D
%
%   See also
%     Patch3D
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
end % end properties


%% Constructor
methods
    function obj = Patch2D(varargin)
    % Constructor for Patch2D class.

        if nargin == 2
            obj.X = varargin{1};
            obj.Y = varargin{2};
        end
    end

end % end constructors


%% Methods specific to Patch2D
methods
    function res = smooth(obj, M)
        % Smooth the patch using local averaging.

        % create convolution vector
        v2 = ones(M, 1) / M;
        
        X2 = padarray(obj.X, [M M], 'replicate');
        X2 = conv2(X2, v2, 'same');
        X2 = X2(M+1:end-M, M+1:end-M);
        
        Y2 = padarray(obj.Y, [M M], 'replicate');
        Y2 = conv2(Y2, v2, 'same');
        Y2 = Y2(M+1:end-M, M+1:end-M);

        % convert result to Patch2D object
        res = Patch2D(X2, Y2, Z2);
    end

    function verts = vertices(obj)
        % Return vertices as a new instance of MultiPoint2D.
        coords = [obj.X(:) obj.Y(:)];
        verts = MultiPoint3D(coords);
    end
    
    function drawSubGrid(varargin)
        % Draw a grid within this patch
        %
        %   drawSubGrid(OBJ, 1) simply displays the boundary of the patch.
        %
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % determines the number of patch tiles
        nTiles = 1;
        if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
            nTiles = varargin{1};
            varargin(1) = [];
        end
        
        % compute indices of curves
        subI = round(linspace(1, size(obj.X,1), nTiles+1));
        subJ = round(linspace(1, size(obj.X,2), nTiles+1));

        % draw 2D curves in first direction
        for i = 1:length(subI)
            x = obj.X(subI(i), :)';
            y = obj.Y(subI(i), :)';

            hx = plot(ax, x, y, varargin{:});
        end
        
        % draw 2D curves in second direction
        for i = 1:length(subJ)
            x = obj.X(:, subJ(i));
            y = obj.Y(:, subJ(i));
            hy = plot(ax, x, y, varargin{:});
        end
        
        if ~isempty(style)
            apply(style, hx);
            apply(style, hy);
        end
    end
end

%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transfo)
        % Apply a geometric transform to this geometry.
        pts = [obj.X(:) obj.Y()];
        pts = transformPoint(transfo, pts);
        dims = size(obj.X);
        res = Patch2D(reshape(pts(:,1), dims), reshape(pts(:,2), dims));
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this patch.
        xmin = min(obj.X(:));
        xmax = max(obj.X(:));
        ymin = min(obj.Y(:));
        ymax = max(obj.Y(:));
        box = Box2D([xmin xmax ymin ymax]);
    end    

    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % extract handle of axis to draw in
        if numel(varargin{1}) == 1 && ishghandle(varargin{1}, 'axes')
            hAx = varargin{1};
            varargin(1) = [];
        else
            hAx = gca;
        end

        % extract the point instance from the list of input arguments
        obj = varargin{1};
        varargin(1) = [];
        
        % extract optional drawing options
        options = {};
        if nargin > 1 && ischar(varargin{1})
            options = varargin;
            varargin = {};
        end
        
        % draw 2D curves in first direction
        for i = 1:size(obj.X, 1)
            x = obj.X(i, :)';
            y = obj.Y(i, :)';

            plot(hAx, x, y, options{:});
        end
        
        % draw 2D curves in second direction
        for j = 1:size(obj.X, 2)
            x = obj.X(:, j);
            y = obj.Y(:, j);
            plot(hAx, x, y, options{:});
        end

        % optionnally add style processing
        if ~isempty(varargin) && isa(varargin{1}, 'Style')
            apply(varargin{1}, hh);
        end
        
        % format output argument
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        res = Patch2D(obj.X * factor, obj.Y * factor);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        res = Patch2D(obj.X + shift(1), obj.Y + shift(2));
    end
    
    function res = rotate(obj, theta)
        % Return a rotated version of this geometry.
        transfo = AffineTransform2D.createRotation(theta);
        res = transform(obj, transfo);
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Patch2D', 'X', obj.X, 'Y', obj.Y);
    end
end
methods (Static)
    function obj = fromStruct(str)
        % Create a new instance from a structure.
        obj = Patch2D(str.X, str.Y);
    end
end

end % end classdef


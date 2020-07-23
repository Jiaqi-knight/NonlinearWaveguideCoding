classdef Box2D < handle
% Bounding box of a planar shape.
%
%   Class Box2D
%   Defined by max extent in each dimension:
%   * XMin, XMax, YMin, YMax.
%
%   Example
%     box = Box2D([5 15 6 14]);
%     figure; axis([0 20 0 20]); hold on
%     draw(box, 'b')
%
%   See also
%     Geometry2D, Box3D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Properties
properties
    XMin;
    XMax;
    YMin;
    YMax;
    
end % end properties


%% Constructor
methods
    function obj = Box2D(varargin)
        % Constructor for Box2D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 1) ~= 1
                error('Creating a box requires an array with one row, not %d', size(var1, 1));
            end
            if size(var1, 2) ~= 4
                error('Creating a box requires an array with four columns, not %d', size(var1, 2));
            end
            data = var1;
        else
            % default box is unit square, with origin as lower-left corner.
            data = [0 1 0 1];
        end
        
        obj.XMin = data(1);
        obj.XMax = data(2);
        obj.YMin = data(3);
        obj.YMax = data(4);
    end

end % end constructors


%% Methods
methods
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        box = Box2D([obj.XMin obj.XMax obj.YMin obj.YMax]);
    end
    
    function varargout = draw(obj, varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % extract style agument if present
        style = [];
        if nargin > 1 && isa(varargin{1}, 'Style')
            style = varargin{1};
            varargin(1) = [];
        end
        
        % draw the box
        tx = [obj.XMin obj.XMax obj.XMax obj.XMin obj.XMin];
        ty = [obj.YMin obj.YMin obj.YMax obj.YMax obj.YMin];
        h = plot(tx, ty, varargin{:});
        
        % eventually apply style
        if ~isempty(style)
            apply(style, h);
        end
        
        % return handle if requested
        if nargout > 0
            varargout = {h};
        end
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of obj geometry.
        factor = varargin{1};
        res = Box2D([obj.XMin obj.XMax obj.YMin obj.YMax] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of obj geometry.
        shift = varargin{1};
        data2 = [obj.XMin obj.XMax obj.YMin obj.YMax] + shift(1, [1 1 2 2]);
        res = Box2D(data2);
    end
    
    function res = rotate(obj, angle) %#ok<STOUT,INUSD>
        % Throw an error as a box can not be rotated.
        error('A box can not be rotated');
    end
end % end methods


%% Serialization methods
methods
    function write(obj, fileName, varargin)
        %WRITE Write box representation into a JSON file.
        % 
        % Requires implementation of the "toStruct" method.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
    
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('type', 'Box2D', ...
            'XMin', obj.XMin, 'XMax', obj.XMax, ...
            'YMin', obj.YMin, 'YMax', obj.YMax);
    end
end

methods (Static)
    function box = read(fileName)
        %READ Read box information from a file in JSON format.
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        box = Box2D.fromStruct(loadjson(fileName));
    end
    
    function box = fromStruct(str)
        % Create a new instance from a structure.
        box = Box2D([str.XMin str.XMax str.YMin str.YMax]);
    end
end

end % end classdef


classdef Sphere3D < Geometry3D
% A sphere, defined by a center and a radius.
%
%   Class Sphere3D
%
%   Example
%   Sphere3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
    % The coordinates of the center.
    Center = [0 0 0];
    
    % The radius of the sphere. Default is 1.
    Radius  = 1;
    
end % end properties


%% Constructor
methods
    function obj = Sphere3D(varargin)
        % Constructor for Sphere3D class
        if nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Sphere3D')
                % copy constructor
                obj.Center = var1.Center;
                obj.Radius  = var1.Radius;
                
            elseif isnumeric(var1)
                % initialize with a 1-by-4 row vector
                obj.Center = var1(1, 1:3);
                obj.Radius  = var1(4);
            end
        elseif nargin == 2
            % center + radius
            
            % initialize center
            var1 = varargin{1};
            if isa(var1, 'Point3D')
                obj.Center = [var1.X var1.Y var1.Z];
            elseif isnumeric(var1)
                obj.Center = var1(1, 1:3);
            else
                error('Can not interpret first argument');
            end
            
            % initialize radius
            obj.Radius = varargin{2};
            
        elseif nargin > 2
            error('Can not initialize sphere');
        end
    end

end % end constructors



%% Methods implementing the Geometry3D interface
methods
    function res = transform(obj, transform) %#ok<STOUT,INUSD>
        % Apply a geometric transform to this geometry.
        error('Method not implemented');
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        
        box = Box3D([...
            obj.Center(1) - obj.Radius obj.Center(1) + obj.Radius ...
            obj.Center(2) - obj.Radius obj.Center(2) + obj.Radius ...
            obj.Center(3) - obj.Radius obj.Center(3) + obj.Radius]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this sphere, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default values for drawing
        nPhi    = 32;
        nTheta  = 16;
        
        % process input options: when a string is found, assumes this is the
        % beginning of options
        options = {'FaceColor', 'g', 'LineStyle', 'none'};
        if length(varargin) == 1
            options = {'FaceColor', varargin{1}, 'LineStyle', 'none'};
        else
            options = [options varargin];
        end
        
        % compute spherical coordinates
        theta   = linspace(0, pi, nTheta+1);
        phi     = linspace(0, 2*pi, nPhi+1);
        
        % convert to cartesian coordinates
        sintheta = sin(theta);
        x = obj.Center(1) + cos(phi') * sintheta * obj.Radius;
        y = obj.Center(2) + sin(phi') * sintheta * obj.Radius;
        z = obj.Center(3) + ones(length(phi),1) * cos(theta) * obj.Radius;
        
        % draw the geometric primitive
        hh = surf(ax, x, y, z, options{:});
        
        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
        
        % format output argument
        if nargout > 0
            h = hh;
        end
    end
end

%% Methods implementing the Geometry3D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        if ~isscalar(factor)
            error('Requires scaling factor to be a scalar');
        end
        center = obj.Center * factor;
        radius = obj.Radius * factor;
        res = Sphere3D(center, radius);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        center = obj.Center + shift;
        res = Sphere3D(center, obj.Radius);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Sphere3D', ...
            'Center', obj.Center, ...
            'Radius', obj.Radius);
    end
end
methods (Static)
    function sphere = fromStruct(str)
        % Create a new instance from a structure.
        sphere = Sphere3D(str.Center, str.Radius);
    end
end

end % end classdef


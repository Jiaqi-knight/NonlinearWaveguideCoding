classdef AffineTransform3D < handle
% A 3D affine transform defined by its matrix.
%
%   Class AffineTransform3D
%
%   Example
%     T = AffineTransform3D.createScaling([3 2 1]);
%
%   See also
%     AffineTransform2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-03,    using Matlab 9.5.0.944444 (R2018b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The coefficients of the transform, as a 1-by-12 numeric array
    % (initialized to identity)
    Coeffs = [1 0 0 0   0 1 0 0   0 0 1 0];
    
end % end properties

%% Static factories
methods (Static)
    function obj = createTranslation(shift, varargin)
        % Create a new 3D affine transform representing a translation.
        % trans = AffineTransform2D.createTranslation([dx dy])
        if isnumeric(shift)
            if all(size(shift) == [1 3])
                dx = shift(1);
                dy = shift(2);
                dz = shift(3);
            elseif isscalar(shift) && isscalar(varargin{1}) && isscalar(varargin{2})
                dx = shift;
                dy = varargin{1};
                dz = varargin{2};
            end
        else
            error('requires numeric input');
        end
         
        obj = AffineTransform3D([1 0 0 dx ; 0 1 0 dy ; 0 0 1 dz; 0 0 0 1]);
    end
    
    function obj = createScaling(factor)
        % Create a new affine transform representing a scaling.
        %
        % Usage:
        % TRANS = AffineTransform2D.createScaling(S);
        % TRANS = AffineTransform2D.createScaling([SX SY SZ]);
        %
        % Example
        %     T = AffineTransform3D.createScaling([3 2 1]);
        
        if isscalar(factor)
            sx = factor;
            sy = factor;
            sz = factor;
        else
            sx = factor(1);
            sy = factor(2);
            sz = factor(3);
        end
        obj = AffineTransform3D([sx 0 0 0   0 sy 0 0   0 0 sz 0]);
    end
    
    function obj = identity()
        % Create a new 3D affine transform representing identity.
        obj = AffineTransform3D([1 0 0 0   0 1 0 0   0 0 1 0]);
    end
end


%% Constructor
methods
    function obj = AffineTransform3D(coeffs)
    % Constructor for AffineTransform3D class.
    
        if nargin < 1
            coeffs = eye(4);
        end
        if isnumeric(coeffs) && all(size(coeffs) == 4)
            coeffs = [coeffs(1,:) coeffs(2,:) coeffs(3,:)];
        end
        if ~(isnumeric(coeffs) && all(size(coeffs) == [1 12]))
            error('requires 1-by-12 or a 4-by-4 numeric array as input');
        end
    
        obj.Coeffs = coeffs;
    end
    
end % end constructors


%% Methods
methods
    function pts2 = transformPoint(obj, pts)
        % Apply transform to a set of coordinates.
        %
        % P2 = transformPoint(T, P)
        % T is the transform object, and P should be a N-by-3 numeric array
        % representing point coordinates.
        % The result P2 has the same size as the input array P.
        %
        
        coeffs = obj.Coeffs;
        pts2 = zeros(size(pts));
        pts2(:,1) = coeffs(1) * pts(:,1) + coeffs(2) * pts(:,2) + coeffs(3) * pts(:,3) + coeffs(4);
        pts2(:,2) = coeffs(5) * pts(:,1) + coeffs(6) * pts(:,2) + coeffs(7) * pts(:,3) + coeffs(8);
        pts2(:,3) = coeffs(9) * pts(:,1) + coeffs(10) * pts(:,2) + coeffs(11) * pts(:,3) + coeffs(12);
    end
    
    function b = isIdentity(obj, varargin)
        % Check if this transfom is equal to identity.
        if isempty(varargin)
            tol = 1e-10;
        else
            tol = varargin{1};
        end
        b = all(abs(obj.Coeffs - [1 0 0 0  0 1 0 0  0 0 1 0]) < tol);
    end
    
    function res = concatenate(obj, obj2)
        % CONCATENATE Concatenate two transforms.
        mat2 = affineMatrix(obj) * affineMatrix(obj2);
        res = AffineTransform3D(mat2);
    end
    
    
    function res = inverse(obj)
        % Return the inverse affine transform.
        res = AffineTransform3D(inv(affineMatrix(obj)));
    end

    function mat = affineMatrix(obj)
        % Get the coefficients of the transform as a 4-by-4 matrix.
        mat = [obj.Coeffs(1:4) ; obj.Coeffs(5:8) ; obj.Coeffs(9:12) ; 0 0 0 1];
    end

end % end methods

%% Overload Matlab computation functions
methods
    function res = mtimes(obj, obj2)
        % Overload the mtime operator.
        mat2 = affineMatrix(obj) * affineMatrix(obj2);
        res = AffineTransform3D(mat2);
    end
end


%% Serialization methods
methods
    function write(obj, fileName, varargin)
        % WRITE Write transform representation into a JSON file.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
    
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct(...
            'Type', 'AffineTransform3D', ...
            'Matrix', [obj.Coeffs(1:4) ; obj.Coeffs(5:8) ; obj.Coeffs(9:12) ; 0 0 0 1]);
    end
end

methods (Static)
    function box = read(fileName)
        %READ Read transform information from a file in JSON format.
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        box = AffineTransform2D.fromStruct(loadjson(fileName));
    end
    
    function transfo = fromStruct(str)
        % Create a new instance from a structure.
        transfo = AffineTransform3D(str.Matrix);
    end
end

end % end classdef


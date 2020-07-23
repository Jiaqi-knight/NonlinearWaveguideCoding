classdef AffineTransform2D < handle
% A 2D affine transform defined by its matrix.
%
%   Class AffineTransform2D
%
%   Example
%     % Create sample geometry and transform, compute the transformed
%     % geometry, and display all results.
%     poly1 = LineString2D([10 10;20 10;20 20;10 20;10 30;20 30]);
%     figure; axis equal; axis([0 50 0 50]);hold on;
%     draw(poly1, 'b');
%     t1 = AffineTransform2D.createTranslation([-10 -20]);
%     t2 = AffineTransform2D.createScaling(2);
%     t3 = AffineTransform2D.createTranslation([25 25]);
%     transfo = t3 * t2 * t1;
%     poly2b = transform(poly1, transfo);
%     draw(poly2b, 'k');
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-03,    using Matlab 9.5.0.944444 (R2018b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % the coefficients of this transform, in order [m00 m01 m02 m10 m11 m12].
    % (initialized to identity)
    Coeffs = [1 0 0  0 1 0];
    
end % end properties


%% Static factories
methods (Static)
    function obj = createTranslation(shift, varargin)
        % Create a new affine transform representing a translation.
        %
        % trans = AffineTransform2D.createTranslation([dx dy])
        if isnumeric(shift)
            if all(size(shift) == [1 2])
                dx = shift(1);
                dy = shift(2);
            elseif isscalar(shift) && isscalar(varargin{1})
                dx = shift;
                dy = varargin{1};
            end
        else
            error('requires numeric input');
        end
         
        obj = AffineTransform2D([1 0 dx 0 1 dy]);
    end
    
    function obj = createRotation(theta)
        % Create a new affine transform representing a rotation.
        %
        % trans = AffineTransform2D.createRotation(angleInRadians)
        cot = cos(theta);
        sit = sin(theta);
        obj = AffineTransform2D([cot -sit 0 sit cot 0]);
    end
    
    function obj = createScaling(factor)
        % Create a new affine transform representing a scaling.
        %
        % Usage:
        % TRANS = AffineTransform2D.createScaling(S);
        % TRANS = AffineTransform2D.createScaling([SX SY]);
        if isscalar(factor)
            sx = factor;
            sy = factor;
        else
            sx = factor(1);
            sy = factor(2);
        end
        obj = AffineTransform2D([sx 0 0   0 sy 0]);
    end
    
    function obj = identity()
        obj = AffineTransform2D([1 0 0   0 1 0]);
    end
end


%% Constructor
methods
    function obj = AffineTransform2D(coeffs)
    % Constructor for AffineTransform2D class.
    
        if nargin < 1
            coeffs = [1 0 0  0 1 0];
        end
        if ~isnumeric(coeffs) 
            error('requires a numeric input');
        end
        if all(size(coeffs) == 3)
            % convert 3x3 matrix to 1x6 row vector (drop last row)
            coeffs = [coeffs(1,1:3) coeffs(2,1:3)];
        end
        if ~all(size(coeffs) == [1 6])
            error('Requires a 1-by-6 numeric array as input');
        end
    
        obj.Coeffs = coeffs;
    end
    
end % end constructors


%% Methods
methods    
    function pts2 = transformPoint(obj, pts)
        % Apply this transform to a set of coordinates.
        %
        % P2 = transformPoint(T, P)
        % T is the transform object, and P should be a N-by-2 numeric array
        % representing point coordinates.
        % The result P2 has the same size as the input array P.
        
        coeffs = obj.Coeffs;
        pts2 = zeros(size(pts));
        pts2(:,1) = coeffs(1) * pts(:,1) + coeffs(2) * pts(:,2) + coeffs(3);
        pts2(:,2) = coeffs(4) * pts(:,1) + coeffs(5) * pts(:,2) + coeffs(6);
    end
    
    function vect2 = transformVector(obj, vect)
        %TRANSFORMVECTOR Apply this affine transform to a vector.
        %
        %   VECT2 = transformVector(T, VECT);
        %   where T is the transformation object, and VECT is the vector to
        %   transform, as a N-by-2 numeric array. 
        %   The result is a N-by-2 numeric array containing the coordinates
        %   of the transformed vector.
        
        % compute new position of vector
        vect2 = zeros(size(vect));
        vect2(:,1) = vect(:,1) * obj.Coeffs(1) + vect(:,2) * obj.Coeffs(2);
        vect2(:,2) = vect(:,1) * obj.Coeffs(4) + vect(:,2) * obj.Coeffs(5);
    end
    
    function b = isIdentity(obj, varargin)
        if isempty(varargin)
            tol = 1e-10;
        else
            tol = varargin{1};
        end
        b = all(abs(obj.Coeffs - [1 0 0 0 1 0]) < tol);
    end
    
    function res = concatenate(obj, obj2)
        %CONCATENATE Concatenate two transforms.
        mat2 = affineMatrix(obj) * affineMatrix(obj2);
        res = AffineTransform2D(mat2);
    end
    
    function res = invert(obj)
        % Return the inverse transform (deprecated).
        % Depreated: should use 'inverse' instead.
        warning('Depreated: should use "inverse" instead.');
        res = AffineTransform2D(inv(affineMatrix(obj)));
    end
    function res = inverse(obj)
        % Return the inverse transform.
        res = AffineTransform2D(inv(affineMatrix(obj)));
    end

    function mat = affineMatrix(obj)
        % Get the coefficients of the transform as a 3-by-3 matrix.
        mat = [obj.Coeffs(1:3) ; obj.Coeffs(4:6) ; 0 0 1];
    end

end % end methods


%% Overload Matlab computation functions
methods
    function res = mtimes(obj, obj2)
        % overload the mtimes method.
        mat2 = affineMatrix(obj) * affineMatrix(obj2);
        res = AffineTransform2D(mat2);
    end
end


%% Serialization methods
methods
    function write(obj, fileName, varargin)
        %WRITE Write transform representation into a JSON file.
        % 
        % Requires implementation of the "toStruct" method.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
    
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct(...
            'Type', 'AffineTransform2D', ...
            'Matrix', [obj.Coeffs(1:3) ; obj.Coeffs(4:6) ; 0 0 1]);
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
        transfo = AffineTransform2D(str.Matrix);
    end
end

end % end classdef


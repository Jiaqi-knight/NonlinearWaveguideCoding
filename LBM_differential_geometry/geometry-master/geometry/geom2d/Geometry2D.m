classdef Geometry2D < Geometry
% Abstract class for planar geometries.
%
%   Class Geometry2D
%
%   Example
%   Geometry2D
%
%   See also
%     Point2D, LineString2D, Polygon2D, MultiPoint2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2014-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.

%% Constructor
methods
    function obj = Geometry2D(varargin)
    % Constructor for Geometry2D class.
    %   (will be called by subclasses)

    end

end % end constructors


%% Abstract Methods
% Declares some methods that will be implemented by subclasses.
methods ( Abstract )
    % Return the bounding box of this geometry.
    box = boundingBox(obj)
    
    %DRAW Draw the current geometry, eventually specifying the style.
    varargout = draw(obj, varargin)

%     % Applies a geometric transform to this geometry
%     res = transform(obj, transform)
end

%% Abstract Methods
% Not sure we will keep these methods...
methods ( Abstract )
    
    % Return a scaled version of this geometry.
    res = scale(obj, varargin)
        
    % Return a translated version of this geometry.       
    res = translate(obj, varargin)
    
    % Return a rotated version of this geometry.
    res = rotate(obj, varargin)
        
end % end methods

end % end classdef


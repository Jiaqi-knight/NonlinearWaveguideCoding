classdef Geometry3D < Geometry
% Abstract class for 3D geometries.
%
%   Class Geometry3D
%
%   Example
%   Geometry3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-02-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRAE - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods
    function obj = Geometry3D(varargin)
    % Constructor for Geometry3D class
    %   (will be called by subclasses)

    end

end % end constructors


%% Abstract Methods
% Declares some methods that will be implemented by subclasses.
methods ( Abstract )
    
    box = boundingBox(obj)
    % Return the 3D bounding box of this geometry.
    
    varargout = draw(obj, varargin)
    %DRAW Draw the current geometry, eventually specifying the style.
end

%% Abstract Methods
% Not sure we will keep these methods...
methods ( Abstract )
    
    res = scale(obj, varargin)
    % Return a scaled version of obj geometry.
        
    res = translate(obj, varargin)
    % Return a translated version of obj geometry.    
    
end % end methods


end % end classdef


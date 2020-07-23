function setupGeometry(varargin)
% Add the different directories of Geometry toolbox to the path.
%
%   Usage:
%   setupGeometry;
%
%   Example
%   setupGeometry;
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-09-16,    using Matlab 9.6.0.1072779 (R2019a)
% Copyright 2011 INRA - Cepia Software Platform.

% extract library path
fileName = mfilename('fullpath');
libDir = fileparts(fileName);

moduleNames = {...
    '.', ...
    'geom2d', ...
    'geom3d'};

disp('Installing Geometry Library');
addpath(libDir);

% add all library modules
for i = 1:length(moduleNames)
    name = moduleNames{i};
    fprintf('Adding module: %-20s', name);
    addpath(fullfile(libDir, name));
    disp(' (ok)');
end


function tests = test_Mesh3D(varargin)
% Test case for the file Mesh3D.
%
%   Test case for the file Mesh3D

%   Example
%   test_Mesh3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-02-07,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRAE - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_createOctahedron(testCase) %#ok<*DEFNU>

mesh = Mesh3D.createOctahedron;
assertTrue(testCase, isa(mesh, 'Mesh3D'));

function test_read_off_mushroom(testCase) %#ok<*DEFNU>
% Test call of function without argument

mesh = Mesh3D.read_off('mushroom.off');
assertTrue(testCase, isa(mesh, 'Mesh3D'));






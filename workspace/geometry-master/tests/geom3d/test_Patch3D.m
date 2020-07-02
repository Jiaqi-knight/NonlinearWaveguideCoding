function tests = test_Patch3D(varargin)
%TEST_PATCH3D  Test case for the file Patch3D
%
%   Test case for the file Patch3D

%   Example
%   test_Patch3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-25,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

x = [1 1 1;2 2 2;3 3 3];
y = [1 2 3;1 2 3;1 2 3];
z = [3 4 3;4 5 4;3 4 3];
patch = Patch3D(x, y, z);
assertTrue(testCase, isa(patch, 'Patch3D'));

function test_Serialize(testCase)
% Test call of function without argument

x = [1 1 1;2 2 2;3 3 3];
y = [1 2 3;1 2 3;1 2 3];
z = [3 4 3;4 5 4;3 4 3];
patch = Patch3D(x, y, z);
str = toStruct(patch);
patch2 = Patch3D.fromStruct(str);
assertEqual(testCase, patch.X, patch2.X);
assertEqual(testCase, patch.Y, patch2.Y);
assertEqual(testCase, patch.Z, patch2.Z);


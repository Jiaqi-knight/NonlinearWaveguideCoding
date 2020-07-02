function tests = test_Point3D(varargin)
%TEST_POINT3D  Test case for the file Point3D
%
%   Test case for the file Point3D

%   Example
%   test_Point3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-07,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_Single(testCase) %#ok<*DEFNU>
% Test call of function without argument

p = Point3D([5 4 3]);
assertEqual(testCase, 5, p.X);
assertEqual(testCase, 4, p.Y);
assertEqual(testCase, 3, p.Z);

function test_Serialize(testCase)
% Test call of function without argument

p = Point3D([5 4 3]);
str = toStruct(p);
p2 = Point3D.fromStruct(str);
assertEqual(testCase, 5, p2.X);
assertEqual(testCase, 4, p2.Y);
assertEqual(testCase, 3, p2.Z);





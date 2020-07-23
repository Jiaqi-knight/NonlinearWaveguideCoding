function tests = test_Sphere3D
%TEST_SPHERE3D  Test case for the file Sphere3D.
%
%   Test case for the file Sphere3D
%
%   Example
%   test_Sphere3D
%
%   See also
%   Sphere3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-12,    using Matlab 9.7.0.1247435 (R2019b) Update 2
% Copyright 2020 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Constructor_TwoArgs(testCase) %#ok<*DEFNU>
% Test call of function without argument

S = Sphere3D([5 4 3], 2);

assertEqual(testCase, S.Center(1), 5);
assertEqual(testCase, S.Center(2), 4);
assertEqual(testCase, S.Center(3), 3);
assertEqual(testCase, S.Radius, 2);


function test_Serialize(testCase)
% Test call of function without argument

S = Sphere3D([5 4 3], 2);
str = toStruct(S);
S2 = Sphere3D.fromStruct(str);
assertEqual(testCase, S.Center(1), S2.Center(1));
assertEqual(testCase, S.Center(2), S2.Center(2));
assertEqual(testCase, S.Center(3), S2.Center(3));
assertEqual(testCase, S.Radius, S2.Radius);

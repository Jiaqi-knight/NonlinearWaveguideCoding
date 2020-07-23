function tests = test_Ellipsoid3D
%TEST_SPHERE3D  Test case for the file Ellipsoid3D.
%
%   Test case for the file Ellipsoid3D
%
%   Example
%   test_Ellipsoid3D
%
%   See also
%   Ellipsoid3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-12,    using Matlab 9.7.0.1247435 (R2019b) Update 2
% Copyright 2020 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Constructor_TwoArgs(testCase) %#ok<*DEFNU>
% Test call of function without argument

E = Ellipsoid3D([30 20 10   5 3 1  30 20 10]);

assertEqual(testCase, E.Center(1), 30);
assertEqual(testCase, E.Center(2), 20);
assertEqual(testCase, E.Center(3), 10);
assertEqual(testCase, E.Radius(1), 5);
assertEqual(testCase, E.Radius(2), 3);
assertEqual(testCase, E.Radius(3), 1);


function test_Serialize(testCase)
% Test call of function without argument

E = Ellipsoid3D([30 20 10   5 3 1  30 20 10]);
str = toStruct(E);
E2 = Ellipsoid3D.fromStruct(str);
assertEqual(testCase, E.Center(1), E2.Center(1));
assertEqual(testCase, E.Center(2), E2.Center(2));
assertEqual(testCase, E.Center(3), E2.Center(3));
assertEqual(testCase, E.Radius(1), E2.Radius(1));
assertEqual(testCase, E.Radius(2), E2.Radius(2));
assertEqual(testCase, E.Radius(3), E2.Radius(3));
assertEqual(testCase, E.EulerAngles(1), E2.EulerAngles(1));
assertEqual(testCase, E.EulerAngles(2), E2.EulerAngles(2));
assertEqual(testCase, E.EulerAngles(3), E2.EulerAngles(3));

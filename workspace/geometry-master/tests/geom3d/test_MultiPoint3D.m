function tests = test_MultiPoint3D(varargin)
%TEST_POINT3D  Test case for the file MultiPoint3D
%
%   Test case for the file MultiPoint3D

%   Example
%   test_MultiPoint3D
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

v = [0 0 0;1 0 0;0 1 0;0 0 1];
pts = MultiPoint3D(v); %#ok<NASGU>



function test_subsref_s1_scalar(testCase) %#ok<*DEFNU>
% test syntax poly(IND) -> return point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = MultiPoint3D(pts);
res = poly(2);
assertTrue(testCase, isa(res, 'Point3D'));

function test_subsref_s1_multi(testCase) %#ok<*DEFNU>
% test syntax poly(INDS) -> return multi-point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = MultiPoint3D(pts);
res = poly([2 3]);
assertTrue(testCase, isa(res, 'MultiPoint3D'));

function test_subsref_s2(testCase) %#ok<*DEFNU>
% test syntax poly(S1, S2) -> return double array

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = MultiPoint3D(pts);
res = poly(1, 1);
assertTrue(testCase, isnumeric(res));
assertTrue(testCase, isscalar(res));

res = poly([2 3], 1);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [2 1]);

res = poly([2 3], [1 2 3]);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [2 3]);

res = poly(:, 1);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [4 1]);

res = poly(2, :);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [1 3]);

res = poly(:, :);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [4 3]);


function test_Serialize(testCase)
% Test call of function without argument

v = [0 0 0;1 0 0;0 1 0;0 0 1];
pts = MultiPoint3D(v);

str = toStruct(pts);

pts2 = MultiPoint3D.fromStruct(str);
assertEqual(testCase, pts.Coords, pts2.Coords);


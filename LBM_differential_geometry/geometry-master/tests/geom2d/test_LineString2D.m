function tests = test_LineString2D
%TEST_LINESTRING2D  Test case for the file LineString2D
%
%   Test case for the file LineString2D
%
%   Example
%   test_LineString2D
%
%   See also
%   LineString2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2020-02-09,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRA - Cepia Software Platform.


tests = functiontests(localfunctions);

function test_Constructor_CoordsArray(testCase) %#ok<*DEFNU>
% Test call of function without argument

pts = [0 0; 10 0;10 10;0 10];
poly = LineString2D(pts);

assertTrue(testCase, isa(poly, 'LineString2D'));


function test_subsref_s1_scalar(testCase) %#ok<*DEFNU>
% test syntax poly(IND) -> return point

pts = [0 0; 10 0;10 10;0 10];
poly = LineString2D(pts);
res = poly(2);
assertTrue(testCase, isa(res, 'Point2D'));

function test_subsref_s1_multi(testCase) %#ok<*DEFNU>
% test syntax poly(INDS) -> return multi-point

pts = [0 0; 10 0;10 10;0 10];
poly = LineString2D(pts);
res = poly([2 3]);
assertTrue(testCase, isa(res, 'MultiPoint2D'));

function test_subsref_s2(testCase) %#ok<*DEFNU>
% test syntax poly(S1, S2) -> return double array

pts = [0 0; 10 0;10 10;0 10];
poly = LineString2D(pts);
res = poly(1, 1);
assertTrue(testCase, isnumeric(res));
assertTrue(testCase, isscalar(res));

res = poly([2 3], 1);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [2 1]);

res = poly([2 3], [1 2]);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [2 2]);

res = poly(:, 1);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [4 1]);

res = poly(2, :);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [1 2]);

res = poly(:, :);
assertTrue(testCase, isnumeric(res));
assertEqual(testCase, size(res), [4 2]);



function test_Length(testCase) %#ok<*DEFNU>
% Test call of function without argument

pts = [0 0; 0 10;10 10;10 0; 20 0; 20 10];
poly = LineString2D(pts);

assertEqual(testCase, length(poly), 50, 'AbsTol', 0.01);


function test_simplify_CircleArc(testCase) %#ok<*DEFNU>

t = linspace(0, pi, 200)';
xt = cos(t);
yt = sin(t);
poly = LineString2D([xt yt]);

poly2 = simplify(poly, .2);

assertTrue(testCase, vertexNumber(poly2) < 10);

assertEqual(testCase, poly2.Coords(1,:), [xt(1) yt(1)], 'AbsTol', 0.01);
assertEqual(testCase, poly2.Coords(end,:), [xt(end) yt(end)], 'AbsTol', 0.01);

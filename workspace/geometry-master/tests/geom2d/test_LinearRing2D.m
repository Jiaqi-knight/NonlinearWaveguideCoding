function tests = test_LinearRing2D
%TEST_LINEARRING2D  Test case for the file LinearRing2D.
%
%   Test case for the file LinearRing2D
%
%   Example
%   test_LinearRing2D
%
%   See also
%   LinearRing2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-02-09,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_CoordsArray(testCase) %#ok<*DEFNU>
% Test call of function without argument

pts = [0 0; 10 0;10 10;0 10];
poly = LinearRing2D(pts);

assertTrue(testCase, isa(poly, 'LinearRing2D'));


function test_subsref_s1_scalar(testCase) %#ok<*DEFNU>
% test syntax poly(IND) -> return point

pts = [0 0; 10 0;10 10;0 10];
poly = LinearRing2D(pts);
res = poly(2);
assertTrue(testCase, isa(res, 'Point2D'));

function test_subsref_s1_multi(testCase) %#ok<*DEFNU>
% test syntax poly(INDS) -> return multi-point

pts = [0 0; 10 0;10 10;0 10];
poly = LinearRing2D(pts);
res = poly([2 3]);
assertTrue(testCase, isa(res, 'MultiPoint2D'));

function test_subsref_s2(testCase) %#ok<*DEFNU>
% test syntax poly(S1, S2) -> return double array

pts = [0 0; 10 0;10 10;0 10];
poly = LinearRing2D(pts);
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

pts = [0 0; 10 0; 10 10; 0 10];
poly = LinearRing2D(pts);

assertEqual(testCase, length(poly), 40, 'AbsTol', 0.01);



function test_simplify_Circle(testCase) %#ok<*DEFNU>

t = linspace(0, 2*pi, 200)';
xt = cos(t);
yt = sin(t);
poly = LinearRing2D([xt yt]);

poly2 = simplify(poly, .2);

assertTrue(testCase, vertexNumber(poly2) < 10);

assertEqual(testCase, poly2.Coords(1,:), [xt(1) yt(1)], 'AbsTol', 0.01);

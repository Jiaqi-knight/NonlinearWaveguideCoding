function tests = test_LineString3D
%TEST_LINESTRING3D  Test case for the file LineString3D
%
%   Test case for the file LineString3D
%
%   Example
%   test_LineString3D
%
%   See also
%   LineString3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-02-13,    using Matlab 9.7.0.1247435 (R2019b) Update 2
% Copyright 2020 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);


function test_subsref_s1_scalar(testCase) %#ok<*DEFNU>
% test syntax poly(IND) -> return point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LineString3D(pts);
res = poly(2);
assertTrue(testCase, isa(res, 'Point3D'));

function test_subsref_s1_multi(testCase) %#ok<*DEFNU>
% test syntax poly(INDS) -> return multi-point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LineString3D(pts);
res = poly([2 3]);
assertTrue(testCase, isa(res, 'MultiPoint3D'));

function test_subsref_s2(testCase) %#ok<*DEFNU>
% test syntax poly(S1, S2) -> return double array

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LineString3D(pts);
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


function test_VertexNumber(testCase) %#ok<*DEFNU>
% Test call of function without argument

ring = LineString3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);
assertEqual(testCase, 4, vertexNumber(ring));


function test_length_square(testCase) %#ok<*DEFNU>
% Test call of function without argument

ring = LineString3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);

perim = length(ring);

exp = 30;
assertEqual(testCase, perim, exp, 'AbsTol', .01);


function test_transform(testCase) 

ring = LineString3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);
trans = AffineTransform3D.createTranslation([3 2 1]);

ring2 = transform(ring, trans);

assertTrue(testCase, isa(ring2, 'LineString3D'));


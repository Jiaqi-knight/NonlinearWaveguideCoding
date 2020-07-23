function tests = test_LinearRing3D(varargin)
% Test case for the file LinearRing3D
%
%   Test case for the file LinearRing3D

%   Example
%   test_Polygon2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-09-11,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);


function test_subsref_s1_scalar(testCase) %#ok<*DEFNU>
% test syntax poly(IND) -> return point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LinearRing3D(pts);
res = poly(2);
assertTrue(testCase, isa(res, 'Point3D'));

function test_subsref_s1_multi(testCase) %#ok<*DEFNU>
% test syntax poly(INDS) -> return multi-point

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LinearRing3D(pts);
res = poly([2 3]);
assertTrue(testCase, isa(res, 'MultiPoint3D'));

function test_subsref_s2(testCase) %#ok<*DEFNU>
% test syntax poly(S1, S2) -> return double array

pts = [0 0 0; 10 0 0;10 10 0;0 10 0];
poly = LinearRing3D(pts);
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

ring = LinearRing3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);
assertEqual(testCase, 4, vertexNumber(ring));


function test_perimeter_square(testCase) %#ok<*DEFNU>
% Test call of function without argument

ring = LinearRing3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);

perim = length(ring);

exp = 40;
assertEqual(testCase, perim, exp, 'AbsTol', .01);


function test_transform(testCase) 

ring = LinearRing3D([0 0 0; 10 0 0; 10 10 0; 0 10 0]);
trans = AffineTransform3D.createTranslation([3 2 1]);

ring2 = transform(ring, trans);

assertTrue(testCase, isa(ring2, 'LinearRing3D'));


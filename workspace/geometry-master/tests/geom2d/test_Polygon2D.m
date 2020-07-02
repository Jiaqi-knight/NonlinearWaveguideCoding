function tests = test_Polygon2D(varargin)
%TEST_POLYGON2D  Test case for the file Polygon2D
%
%   Test case for the file Polygon2D

%   Example
%   test_Polygon2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2018-09-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor(testCase) %#ok<*DEFNU>
% Test call of function without argument

vertices = [10 10; 20 10; 20 20; 10 20];
poly = Polygon2D(vertices);
assertEqual(testCase, 4, size(poly.Coords, 1));

function test_perimeter_square(testCase) %#ok<*DEFNU>
% Test call of function without argument

vertices = [10 10; 20 10; 20 20; 10 20];
poly = Polygon2D(vertices);
exp = 40;

perim = perimeter(poly);

assertEqual(testCase, perim, exp, 'AbsTol', .01);

function test_area_square(testCase) %#ok<*DEFNU>
% Test call of function without argument

vertices = [10 10; 20 10; 20 20; 10 20];
poly = Polygon2D(vertices);
exp = 100;

a = area(poly);

assertEqual(testCase, a, exp, 'AbsTol', .01);

function test_transform(testCase) 

vertices = [10 10; 20 10; 20 20; 10 20];
poly = Polygon2D(vertices);
trans = AffineTransform2D.createRotation(pi/3);
poly2 = transform(poly, trans);

assertTrue(testCase, isa(poly2, 'Polygon2D'));


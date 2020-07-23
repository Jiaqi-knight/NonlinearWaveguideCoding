function tests = test_Point2D
%TEST_POINT2D  Test case for the file Point2D
%
%   Test case for the file Point2D
%
%   Example
%   test_Point2D
%
%   See also
%   Point2D
%
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-03,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_Single(testCase) %#ok<*DEFNU>
% Test call of function without argument

p = Point2D([3 4]);
assertEqual(testCase, 3, p.X);
assertEqual(testCase, 4, p.Y);


function test_subsref(testCase)

p = Point2D([3 4]);
assertEqual(testCase, 3, p(1));
assertEqual(testCase, 4, p(2));


function test_Serialize(testCase)
% Test call of function without argument

p = Point2D([3 4]);
str = toStruct(p);
p2 = Point2D.fromStruct(str);
assertEqual(testCase, 3, p2.X);
assertEqual(testCase, 4, p2.Y);


function test_draw_simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

pt = Point2D([20 10]); 

f = figure;
h = draw(pt);
assertTrue(testCase, ishandle(h));
close(f);


function test_draw_drawOptions(testCase) %#ok<*DEFNU>
% Test call of function without argument

pt = Point2D([20 10]); 

f = figure;
h = draw(pt, 'Color', 'b', 'Marker', '+');
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.Marker, '+');
close(f);

function test_draw_style(testCase) %#ok<*DEFNU>
% Test call of function without argument

pt = Point2D([20 10]); 
style = Style('MarkerVisible', true, 'MarkerColor', 'b', 'MarkerStyle', '+');

f = figure;
h = draw(pt, style);

assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.Marker, '+');

close(f);

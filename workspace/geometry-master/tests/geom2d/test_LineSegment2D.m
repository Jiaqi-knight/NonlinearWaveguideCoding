function tests = test_LineSegment2D
%TEST_LINESEGMENT2D  Test case for the file LineSegment2D
%
%   Test case for the file LineSegment2D
%
%   Example
%   test_LineSegment2D
%
%   See also
%   LineSegment2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_TwoPoints(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);

seg = LineSegment2D(p1, p2);

assertTrue(testCase, isa(seg, 'LineSegment2D'));

function test_Constructor_Empty(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg = LineSegment2D();

assertTrue(testCase, isa(seg, 'LineSegment2D'));

function test_Constructor_Copy(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);

seg1 = LineSegment2D(p1, p2);
seg2 = LineSegment2D(seg1);

assertTrue(testCase, isa(seg2, 'LineSegment2D'));


function test_distancePoint(testCase) %#ok<*DEFNU>

% create line segment with slope 3/4
p1 = Point2D(10, 10);
p2 = Point2D(50, 40);
seg = LineSegment2D(p1, p2);

% points around P1
pt1 = Point2D(10+3, 10-4);
assertEqual(testCase, seg.distance(pt1), 5.0, 'AbsTol', 0.01);
pt2 = Point2D(10-4, 10-3);
assertEqual(testCase, seg.distance(pt2), 5.0, 'AbsTol', 0.01);
pt3 = Point2D(10-3, 10+4);
assertEqual(testCase, seg.distance(pt3), 5.0, 'AbsTol', 0.01);

% points around P1
pt4 = Point2D(50+3, 40-4);
assertEqual(testCase, seg.distance(pt4), 5.0, 'AbsTol', 0.01);
pt5 = Point2D(50+4, 40+3);
assertEqual(testCase, seg.distance(pt5), 5.0, 'AbsTol', 0.01);
pt6 = Point2D(50-3, 40+4);
assertEqual(testCase, seg.distance(pt6), 5.0, 'AbsTol', 0.01);


function test_draw_simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);
seg = LineSegment2D(p1, p2);

f = figure;
h = draw(seg);
assertTrue(testCase, ishandle(h));
close(f);


function test_draw_drawOptions(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);
seg = LineSegment2D(p1, p2);

f = figure;
h = draw(seg, 'Color', 'b', 'LineWidth', 2);
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);
close(f);

function test_draw_axis(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);
seg = LineSegment2D(p1, p2);

f = figure;
ax = gca;

h = draw(ax, seg, 'Color', 'b', 'LineWidth', 2);
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);
close(f);

function test_draw_style(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);
seg = LineSegment2D(p1, p2);
style = Style('LineColor', 'b', 'LineWidth', 2);

f = figure;
h = draw(seg, style);

assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);

close(f);

function tests = test_MultiCurve2D
%TEST_MULTICURVE2D  Test case for the file MultiCurve2D
%
%   Test case for the file MultiCurve2D
%
%   Example
%   test_MultiCurve2D
%
%   See also
%   MultiCurve2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_constructor_init(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg1 = LineSegment2D(Point2D(4, 5), Point2D(6,5));
seg2 = LineSegment2D(Point2D(5, 4), Point2D(5, 6));
circ = Circle2D([7, 3], 2);

mc = MultiCurve2D(seg1, seg2, circ);

assertTrue(testCase, isa(mc, 'MultiCurve2D'));

function test_draw(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg1 = LineSegment2D(Point2D(4, 5), Point2D(6,5));
seg2 = LineSegment2D(Point2D(5, 4), Point2D(5, 6));
circ = Circle2D([7, 3], 2);
mc = MultiCurve2D(seg1, seg2, circ);

fig = figure; axis equal; axis([0 10 0 10]); hold on;
draw(mc, 'b');
close(fig);


function test_Serialize(testCase)
% Test call of function without argument

seg1 = LineSegment2D(Point2D(4, 5), Point2D(6,5));
seg2 = LineSegment2D(Point2D(5, 4), Point2D(5, 6));
circ = Circle2D([7, 3], 2);
mc = MultiCurve2D(seg1, seg2, circ);

str = toStruct(mc);
mc2 = Geometry.fromStruct(str);

assertTrue(isa(mc2, 'MultiCurve2D'));
assertEqual(testCase, length(mc.Curves), length(mc2.Curves));


function test_toStruct(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg1 = LineSegment2D(Point2D(4, 5), Point2D(6,5));
seg2 = LineSegment2D(Point2D(5, 4), Point2D(5, 6));
circ = Circle2D([7, 3], 2);
mc = MultiCurve2D(seg1, seg2, circ);

str = toStruct(mc);

assertEqual(testCase, length(str), 1);

function tests = test_StraightLine2D
%TEST_StraightLine2D  Test case for the file StraightLine2D
%
%   Test case for the file StraightLine2D
%
%   Example
%   test_StraightLine2D
%
%   See also
%   StraightLine2D
%
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-03,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_Single(testCase) %#ok<*DEFNU>
% Test call of function without argument

p = StraightLine2D([3 4], [6 5]);
assertEqual(testCase, [3 4], p.Origin);
assertEqual(testCase, [3 1], p.Direction);


function test_clip(testCase)

box = Box2D([0 10 0 10]);
line = StraightLine2D([5 5 1 .1]);
res = clip(line, box);
assertTrue(testCase, isa(res, 'LineSegment2D'));


function test_Serialize(testCase)
% Test call of function without argument

line = StraightLine2D([3 4], [6 5]);
str = toStruct(line);
line2 = StraightLine2D.fromStruct(str);
assertEqual(testCase, line.Origin, line2.Origin);
assertEqual(testCase, line.Direction, line2.Direction);


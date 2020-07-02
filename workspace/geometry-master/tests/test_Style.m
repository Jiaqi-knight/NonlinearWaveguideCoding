function tests = test_Style(varargin)
%TEST_STYLE  Test case for the file Style
%
%   Test case for the file Style

%   Example
%   test_Style
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2018-09-19,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_EmptyConstructor(testCase) %#ok<*DEFNU>
% Test call of function without argument

style = Style();
style.LineWidth;


function test_ArgumentConstructor(testCase) %#ok<*DEFNU>
% Test call of function without argument

style = Style('LineStyle', '-', 'linewidth', 2);
style.LineWidth;


function test_read(testCase)
% Test call of function without argument

style = Style.read('style1.style');
assertTrue(testCase, isa(style, 'Style')); 
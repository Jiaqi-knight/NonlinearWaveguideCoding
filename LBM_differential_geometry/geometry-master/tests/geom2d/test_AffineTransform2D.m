function tests = test_AffineTransform2D(varargin)
%TEST_AFFINETRANSFORM2D  Test case for the file AffineTransform2D
%
%   Test case for the file AffineTransform2D

%   Example
%   test_AffineTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-20,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_EmptyConstructor(testCase) %#ok<*DEFNU>
% Test call of function without argument
trans = AffineTransform2D();
assertTrue(testCase, isa(trans, 'AffineTransform2D'));

function test_createTranslation(testCase) %#ok<*DEFNU>

trans = AffineTransform2D.createTranslation([20 10]);
assertTrue(testCase, isa(trans, 'AffineTransform2D'));


function test_createRotation(testCase)
trans = AffineTransform2D.createRotation(pi/3);
assertTrue(testCase, isa(trans, 'AffineTransform2D'));


function test_isIdentity_true(testCase) 
trans = AffineTransform2D();
assertTrue(testCase, isIdentity(trans));

function test_isIdentity_false(testCase) 
trans = AffineTransform2D.createRotation(pi/3);
assertFalse(testCase, isIdentity(trans));


function test_toAndFromStruct(testCase)
trans = AffineTransform2D.createTranslation([30 20]);
str = toStruct(trans);
trans2 = AffineTransform2D.fromStruct(str);
assertEqual(testCase, trans.Coeffs, trans2.Coeffs);

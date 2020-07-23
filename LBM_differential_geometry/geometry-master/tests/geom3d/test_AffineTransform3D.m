function tests = test_AffineTransform3D(varargin)
%TEST_AFFINETRANSFORM3D  Test case for the file AffineTransform3D
%
%   Test case for the file AffineTransform3D

%   Example
%   test_AffineTransform3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-05-26,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_EmptyConstructor(testCase) %#ok<*DEFNU>
% Test call of function without argument
trans = AffineTransform3D();
assertTrue(testCase, isa(trans, 'AffineTransform3D'));

function test_createTranslation(testCase) %#ok<*DEFNU>

trans = AffineTransform3D.createTranslation([30 20 10]);
assertTrue(testCase, isa(trans, 'AffineTransform3D'));


% function test_createRotation(testCase)
% trans = AffineTransform3D.createRotation(pi/3);
% assertTrue(testCase, isa(trans, 'AffineTransform3D'));


function test_isIdentity_true(testCase) 
trans = AffineTransform3D();
assertTrue(testCase, isIdentity(trans));

function test_isIdentity_false(testCase) 
trans = AffineTransform3D.createTranslation([5 4 3]);
assertFalse(testCase, isIdentity(trans));


function test_toAndFromStruct(testCase)
trans = AffineTransform3D.createTranslation([5 4 3]);
str = toStruct(trans);
trans2 = AffineTransform3D.fromStruct(str);
assertEqual(testCase, trans.Coeffs, trans2.Coeffs);

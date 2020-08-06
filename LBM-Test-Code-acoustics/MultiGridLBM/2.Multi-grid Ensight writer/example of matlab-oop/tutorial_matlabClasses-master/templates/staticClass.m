classdef staticClass
    % Filename: staticClass.m
    % Author:   Samuel Acuña
    % Date:     19 Dec 2018
    % Description:
    % Sample static class of functions.
    % 
    %
    % Example Usage:
    %       staticClass.double(5)       % = 10
    %       staticClass.squared(5)      % = 25
    %       staticClass.doubleSquare(5) % = 100
    %       staticClass.addConstant(5)  % = 505
    
    properties (Constant)
        % properties are optional, but often helpful
        prop1 = 500;
    end
    
    methods (Static) % collection of your assorted functions. 
        % Don't need a constructor function.
        function y = double(x)
            % function info here
            y = 2*x;
        end
        function y = squared(x)
            % function info here
            y = x^2;
        end
        function z = doubleSquare(x)
            % function info here
            y = staticClass.double(x);
            z = staticClass.squared(y);
        end
        function y = addConstant(x)
            y = x + staticClass.prop1;
        end
    end %methods
end %classdef
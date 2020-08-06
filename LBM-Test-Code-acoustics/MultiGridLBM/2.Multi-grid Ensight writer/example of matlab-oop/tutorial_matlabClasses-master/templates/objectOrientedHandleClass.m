classdef objectOrientedHandleClass < handle
    % Filename: objectOrientedHandleClass.m
    % Author:   Samuel Acuña
    % Date:     19 Dec 2018
    % Description:
    % Sample objected oriented class, implementing a handle type
    % 
    %
    % Example Usage:
    %       x = objectOrientedHandleClass(5) 
    %       x.double()       % = 10
    %       x.squared()      % = 100
    %       x.doubleSquare() % = 40000
    %       x.addValue(500)  % = 40500
    
    %%%%%%%%%%%%%
    % PUBLIC PROPERTIES
    properties
        % properties are defined in the constructor function
        prop1; 
    end
    
    methods 
        %%%%%%%%%%%%
        % CONSTRUCTOR
        % constructor function must have same name as class. This is the
        % only function that needs to explicitly output the object.
        function obj = objectOrientedHandleClass(input)
            obj.prop1 = input;
        end
        
        %%%%%%%%%%%%
        % SPECIAL OBJECT METHODS
        function double(obj)
            % function info here
            obj.prop1 = 2*obj.prop1;
        end
        function squared(obj)
            % function info here
             obj.prop1 =  (obj.prop1)^2;
        end
        function doubleSquare(obj)
            % function info here
            obj.double();
            obj.squared();
        end
        function addValue(obj,input)
            obj.prop1 = obj.prop1 + input;
        end
        
        %%%%%%%%%%%%
        % SET FUNCTIONS (not always necessary, but can be helpful)
        
        %%%%%%%%%%%%
        % GET FUNCTIONS (necessary for dependent properties)
        
    end %methods
end %classdef
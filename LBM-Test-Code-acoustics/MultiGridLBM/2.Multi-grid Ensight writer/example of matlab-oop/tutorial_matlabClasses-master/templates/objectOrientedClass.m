classdef objectOrientedClass
    % Filename: objectOrientedClass.m
    % Author:   Samuel Acuña
    % Date:     19 Dec 2018
    % Description:
    % Sample objected oriented class
    % 
    %
    % Example Usage:
    %       x = objectOrientedClass(5) 
    %       x = x.double()       % = 10
    %       x = x.squared()      % = 100
    %       x = x.doubleSquare() % = 40000
    %       x = x.addValue(500)  % = 40500
    
    %%%%%%%%%%%%%
    % PUBLIC PROPERTIES
    properties
        % properties are defined in the constructor function
        prop1; 
    end

    
    methods 
        %%%%%%%%%%%%
        % CONSTRUCTOR
        % constructor function must have same name as class
        function obj = objectOrientedClass(input)
            obj.prop1 = input;
        end
        
        %%%%%%%%%%%%
        % SPECIAL OBJECT METHODS
        function obj = double(obj)
            % function info here
            obj.prop1 = 2*obj.prop1;
        end
        function obj = squared(obj)
            % function info here
             obj.prop1 =  (obj.prop1)^2;
        end
        function obj = doubleSquare(obj)
            % function info here
            obj = obj.double();
            obj = obj.squared();
        end
        function obj = addValue(obj,input)
            obj.prop1 = obj.prop1 + input;
        end
        
        %%%%%%%%%%%%
        % SET FUNCTIONS (not always necessary, but can be helpful)
        
        %%%%%%%%%%%%
        % GET FUNCTIONS (necessary for dependent properties)
        
    end %methods
end %classdef

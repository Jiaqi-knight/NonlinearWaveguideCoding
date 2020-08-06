classdef weapon
    % Filename: weapon.m
    % Author:   Samuel Acuña
    % Date:     18 Dec 2018
    %
    % Project idea: create a video game
    % Premise: knights fight each other
    %
    % Example 4: matlab objected oriented programming
    %

    
    %%%%%%%%%%%%
    % PUBLIC PROPERTIES
    properties
        type; % type of weapon
        base_attack; % amount of damage done when used
    end
    
    methods 
        function obj = weapon(typeNum) % MUST HAVE SAME NAME AS CLASS
            % CONSTRUCTOR FUNCTION. Goal: instantiate the undefined
            % properties. If user doesnt specify when constructing this
            % object, provide default values.
            % set the armour type from a discrete set of options
            if nargin == 0
                typeNum = randi(5);
            end
            
            switch typeNum
                case 1
                    obj.type = 'sword';
                    obj.base_attack = 5;
                case 2
                    obj.type = 'spear';
                    obj.base_attack = 6;
                case 3
                    obj.type = 'mace';
                    obj.base_attack = 7;
                case 4
                    obj.type = 'lance';
                    obj.base_attack = 8;
                case 5
                    obj.type = 'rubber chicken';
                    obj.base_attack = 1;
                otherwise
                    error('Weapon type number must be between 1 and 5')
            end
        end
    end %methods
end
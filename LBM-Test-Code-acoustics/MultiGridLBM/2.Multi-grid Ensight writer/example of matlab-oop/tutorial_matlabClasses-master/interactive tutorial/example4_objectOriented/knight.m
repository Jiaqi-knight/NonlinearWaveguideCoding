classdef knight
    % Filename: knight.m
    % Author:   Samuel Acuña
    % Date:     18 Dec 2018
    %
    % Project idea: create a video game
    % Premise: knights fight each other
    %
    % Example 4: matlab objected oriented programming
    %
    % for more on property attributes, see:
    % https://www.mathworks.com/help/matlab/matlab_oop/property-attributes.html
    %
    % for more on method attributes, see:
    % https://www.mathworks.com/help/matlab/matlab_oop/method-attributes.html
    %
    
    
    %%%%%%%%%%%%%
    % PUBLIC PROPERTIES
    properties (GetAccess = 'public', SetAccess = 'private')
        % public read access, but private write access.
        name; % name of the knight
        health; % current health
        base_health; % starting level of health
        attack; % current attack power of this knight
        weapon; % weapon equipped
        nHitsRecieved; % how many hits has this knight suffered?
    end
    properties (Dependent) % these properties are calculated from other properties
        isDead; % has this knight died? (depends on health)
    end
    %%%%%%%%%%%%
    % HIDDEN PROPERTIES
    properties (Hidden) % hidden properties
        % gamplay values:
        bonus; % percent bonus in attack after flexing or roaring
        penalty; % percent penalty in attack after showing Off or talking smack
        criticalHit_modifier; % how much a critical hit increases damage taken
        odds_of_criticalHit; % percent chance of taking a critical hit when attacked
    end
    
    
    %%%%%%%%%%%%%
    % PUBLIC METHODS
    methods
        %%%%%%%%%%%%
        % CONSTRUCTOR
        function obj = knight(name)
            % CONSTRUCTOR FUNCTION. Goal: instantiate the undefined
            % properties. If user doesnt specify when constructing this
            % object, provide default values.
            
            obj.name = name; % note: (see set.name) checks that 'name' is a string
            
            obj.base_health = randi(30)+70; % random starting health
            obj.health = obj.base_health;
            obj.weapon = weapon();
            obj.attack = obj.weapon.base_attack;
            obj.nHitsRecieved = 0;
            obj.bonus = 0.1*randi(3); % random variation in bonus stat
            obj.penalty = 0.1*randi(2); % random variation in penalty stat
            obj.criticalHit_modifier = 1.2+0.1*randi(7); % random variation in criticalHit_modifier
            obj.odds_of_criticalHit = 0.10+0.01*randi(40); % random variation in odds of critical hit
            
            disp(['You created a knight named ' obj.name '.']);
        end
        
        %%%%%%%%%%%%
        % SPECIAL OBJECT METHODS
        function obj = getsHit(obj,knight2)
            damage = knight2.attack; % pull damage values from other knight's attack
            disp([obj.name ' was attacked.']);
            damage = obj.isCriticalHit(damage); % see if critical hit, if so increase damage
            disp(['   Recieved ' num2str(damage) ' damage.']);
            obj.health = obj.health - damage; % damage reduces health
            obj.nHitsRecieved = obj.nHitsRecieved + 1; % increment hit counter
            if obj.isDead
                obj.health = 0;
                disp([obj.name ' is dead.']);
                disp('Game Over.')
            end
        end
        function obj = flexes(obj)
            disp([obj.name ' flexes, and their attack increase by ' num2str(obj.bonus*100) '%']);
            obj.attack = (1+obj.bonus)*obj.attack;
        end
        function obj = roars(obj)
            disp([obj.name ' roars, and their attack increase by ' num2str(obj.bonus*100) '%']);
            obj.attack = (1+obj.bonus)*obj.attack;
        end
        function obj = showsOff(obj)
            disp([obj.name ' is showing off, and their attack decreases by ' num2str(obj.penalty*100) '%']);
            obj.attack = (1-obj.penalty)*obj.attack;
        end
        function obj = talksSmack(obj)
            disp([obj.name ' is talking smack, and their attack decreases by ' num2str(obj.penalty*100) '%']);
            obj.attack = (1-obj.penalty)*obj.attack;
        end
        
        %%%%%%%%%%%%
        % SET FUNCTIONS (not always necessary, but can be helpful)
        function obj = set.name(obj,name)
            % sets name of knight, but checks it is a string first
            if ~ischar(name)
                error(['Knight''s name must be a string.']);
            end
            obj.name = name; % set name
        end
        
        %%%%%%%%%%%%
        % GET FUNCTIONS (necessary for dependent properties)
        function TF = get.isDead(obj)
            TF = 0; % default is alive
            if obj.health <= 0
                TF = 1;
            end
        end
    end %methods
    
    %%%%%%%%%%%%%
    % PRIVATE METHODS
    methods (Access = private)
        function newDamage = isCriticalHit(obj,damage)
            % add chance of taking a critical hit
            if rand(1) < obj.odds_of_criticalHit
                disp('It was a critical hit!');
                newDamage = obj.criticalHit_modifier*damage;
            else
                newDamage = damage;
            end
        end
    end %methods
end %classdef
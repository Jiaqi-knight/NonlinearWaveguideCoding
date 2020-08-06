classdef knight
    % Filename: knight.m
    % Author:   Samuel Acuña
    % Date:     18 Dec 2018
    %
    % Project idea: create a video game
    % Premise: knights fight each other
    %
    % Example 3: matlab class with static methods and properties
    
    
    %%%%%%%%%%%%%
    % PROPERTIES
    properties (Constant)
        bonus = 0.20;
        penalty = 0.10;
        criticalHit_modifier = 1.5;
        odds_of_criticalHit = 0.2;
    end
    
    %%%%%%%%%%%%%
    % PUBLIC METHODS
    methods (Static)
        function k = getsHit(k,damage)
            disp([k.name ' was attacked.']);
            damage = knight.isCriticalHit(damage); % see if critical hit, if so increase damage
            disp(['   Recieved ' num2str(damage) ' damage.']);
            k.health = k.health - damage; % damage reduces health
            if knight.isDead(k)
                disp('Game Over.')
            end
        end
        function k = flexes(k)
            disp([k.name ' flexes, and their attack increase by ' num2str(knight.bonus*100) '%']);
            k.attack = (1+knight.bonus)*k.attack;
        end
        function k = roars(k)
            disp([k.name ' roars, and their attack increase by ' num2str(knight.bonus*100) '%']);
            k.attack = (1+knight.bonus)*k.attack;
        end
        function k = showsOff(k)
            disp([k.name ' is showing off, and their attack decreases by ' num2str(knight.penalty*100) '%']);
            k.attack = (1-knight.penalty)*k.attack;
        end
        function k = talksSmack(k)
            disp([k.name ' is talking smack, and their attack decreases by ' num2str(knight.penalty*100) '%']);
            k.attack = (1-knight.penalty)*k.attack;
        end
    end %methods
   
    %%%%%%%%%%%%%
    % PRIVATE METHODS
    methods (Static, Access = private)
        function newDamage = isCriticalHit(damage)
            % add chance of getting a critical hit
            if rand(1) < knight.odds_of_criticalHit
                disp('It was a critical hit!');
                newDamage = knight.criticalHit_modifier*damage;
            else
                newDamage = damage;
            end
        end
        function TF = isDead(k) % check if knight is dead
            TF = 0; % default is alive
            if k.health <= 0
                disp([k.name ' is dead.']);
                TF = 1;
            end
        end
    end %methods
end %classdef
classdef knight
    % Filename: knight.m
    % Author:   Samuel Acuña
    % Date:     18 Dec 2018
    %
    % Project idea: create a video game
    % Premise: knights fight each other
    %
    % Example 2: matlab class with static methods
    %
    
    %%%%%%%%%%%%%
    % PUBLIC METHODS
    methods (Static)
        function k = getsHit(k,damage)
            disp([k.name ' was attacked and recieved ' num2str(damage) ' damage']);
            k.health = k.health - damage; % damage reduces health
            if knight.isDead(k)
                disp('Game Over.')
            end
        end
        function k = flexes(k)
            bonus = 0.20;
            disp([k.name ' flexes, and their attack increase by ' num2str(bonus*100) '%']);
            k.attack = (1+bonus)*k.attack;
        end
        function k = roars(k)
            bonus = 0.20;
            disp([k.name ' roars, and their attack increase by ' num2str(bonus*100) '%']);
            k.attack = (1+bonus)*k.attack;
        end
        function k = showsOff(k)
            penalty = 0.10;
            disp([k.name ' is showing off, and their attack decreases by ' num2str(penalty*100) '%']);
            k.attack = (1-penalty)*k.attack;
        end
        function k = talksSmack(k)
            penalty = 0.10;
            disp([k.name ' is talking smack, and their attack decreases by ' num2str(penalty*100) '%']);
            k.attack = (1-penalty)*k.attack;
        end
    end %methods
    %%%%%%%%%%%%%
    % PRIVATE METHODS
    methods (Static, Access = private)
        function TF = isDead(k) % check if knight is dead
            TF = 0; % default is alive
            if k.health <= 0
                disp([k.name ' is dead.']);
                TF = 1;
            end
        end
    end %methods
end %classdef
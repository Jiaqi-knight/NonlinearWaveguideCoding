function knight = flexes(knight)
% Filename: flexes.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
%
% Project idea: create a video game
% Premise: knights fight each other
%
% Example 2: matlab class with static methods
%

bonus = 0.20;
disp([knight.name ' flexes, and their attack increase by ' num2str(bonus*100) '%']);
knight.attack = (1+bonus)*knight.attack;

end
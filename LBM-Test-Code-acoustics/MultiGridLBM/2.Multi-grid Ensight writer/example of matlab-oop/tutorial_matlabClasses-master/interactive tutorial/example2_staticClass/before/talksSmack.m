function knight = talksSmack(knight)
% Filename: showsOff.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
%
% Project idea: create a video game
% Premise: knights fight each other
%
% Example 2: matlab class with static methods

penalty = 0.10;
disp([knight.name ' is talking smack, and their attack decreases by ' num2str(penalty*100) '%']);
knight.attack = (1-penalty)*knight.attack;


end
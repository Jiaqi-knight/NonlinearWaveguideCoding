function TF = isDead(knight) % check if knight is dead
% Filename: isDead.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
%
% Project idea: create a video game
% Premise: knights fight each other
%
% Example 2: matlab class with static methods
%
% Description:
% checks if knight is dead

TF = 0; % default is alive
if knight.health <= 0
    disp([knight.name ' is dead.']);
    TF = 1;
end

end
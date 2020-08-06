function knight = getsHit(knight,damage)
% Filename: getsHit.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
%
% Project idea: create a video game
% Premise: knights fight each other
%
% Example 2: matlab class with static methods
%

 disp([knight.name ' was attacked and recieved ' num2str(damage) ' damage']);
 knight.health = knight.health - damage; % damage reduces health
 
 if isDead(knight)
     disp('Game Over.')
 end
 
end
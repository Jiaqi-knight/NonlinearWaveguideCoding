% Filename: knightGame.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
% 
% Project idea: create a video game
% Premise: knights fight each other
%
% Example 2: matlab class with static methods
%
%
% NOTE: this code works best by running individual 
% sections (or cells) individually. This can be done by 
% pressing "run section" or "run and advance" in the 
% editor tab (i recommend memorizing the keyboard shortcut).
% for more information, see: https://www.mathworks.com/help/matlab/matlab_prog/run-sections-of-programs.html

clear; close all; clc;

%%%%%%%%%%%%%%
%% Create the characters.
k1.name = 'Simon';
k1.health = 50;
k1.attack = 5;
k1.weapon = 'sword';

k2.name = 'Suzy';
k2.health = 60;
k2.attack = 7;
k2.weapon = 'spear';
% return

%%%%%%%%%%%%%%
%% show knights stats
k1
k2

%%%%%%%%%%%%%%
%% simulate fight actions
k1 = knight.getsHit(k1,5); % what happens if you remove the "k1 = " from the front?
%%
k1 = knight.getsHit(k1,k2.attack);
%%
k2 = knight.getsHit(k2,k1.attack);
%%
k1 = knight.flexes(k1);
%%
k1 = knight.roars(k1);
%%
k2 = knight.showsOff(k2);
%%
k1 = knight.talksSmack(k1);
%%
isDead(k1)


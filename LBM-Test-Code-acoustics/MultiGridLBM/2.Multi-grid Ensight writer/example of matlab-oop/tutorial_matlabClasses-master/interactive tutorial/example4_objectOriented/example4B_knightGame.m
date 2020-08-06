% Filename: example4B_knightGame.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
% 
% Project idea: create a video game
% Premise: knights fight each other to the death
%
% Example 4: matlab objected oriented programming
%
% NOTE: this code works best by running individual 
% sections (or cells) individually. This can be done by 
% pressing "run section" or "run and advance" in the 
% editor tab (i recommend memorizing the keyboard shortcut).
% for more information, see: https://www.mathworks.com/help/matlab/matlab_prog/run-sections-of-programs.html

clear; close all; clc;
% return

%%%%%%%%%%%%%%
%% explore the knight class
k1 = knight('Simon') %creates knight with random stats
k2 = knight('Suzy') 

%% try to given knight invalid name
% k1 = knight(555)

%% look at object of objects
clc;
k1
k1.weapon
k2
k2.weapon


%%%%%%%%%%%%%%
%% show knights stats
k1
k2

%%%%%%%%%%%%%%
%% simulate fight actions
k1 = k1.getsHit(k2); % what happens if you remove the "k = "?
%%
k2 = k2.getsHit(k1);
%%
k1 = k1.flexes();
%%
k1 = k1.roars();
%%
k2 = k2.showsOff();
%%
k1 = k1.talksSmack();



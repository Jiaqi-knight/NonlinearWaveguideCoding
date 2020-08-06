% Filename: knightGame.m
% Author:   Samuel Acuña
% Date:     18 Dec 2018
% 
% Project idea: create a video game
% Premise: knights fight each other to the death
%
% Example 5: matlab objected oriented programming with handles
%
% Handles only explicitly return an object for the constructor function.
% Any changes to the object are saved automatically (dont have to save
% over another copy)
%
% for more on implementing a handle class, see
% https://www.mathworks.com/help/matlab/ref/handle-class.html
%
% NOTE: this code works best by running individual 
% sections (or cells) individually. This can be done by 
% pressing "run section" or "run and advance" in the 
% editor tab (i recommend memorizing the keyboard shortcut).
% for more information, see: https://www.mathworks.com/help/matlab/matlab_prog/run-sections-of-programs.html

clear; close all; clc;
% return

%%%%%%%%%%%%%%
%% set up knights
k1 = knight('Simon'); %creates knight with random stats
k2 = knight('Suzy'); 

%%%%%%%%%%%%%%
%% show knights stats
k1
k2

%%%%%%%%%%%%%%
%% simulate fight actions
% k1 = k1.getsHit(k2); % what happens if you DON'T remove the "k = "?
%%
k2.getsHit(k1);
%%
k1.flexes();
%%
k1.roars();
%%
k2.showsOff();
%%
k1.talksSmack();



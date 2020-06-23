%spherium
% Top level function which adds the spherium directory to the path with
% subolders, and runs the top level GUI spherium_gui.m

function spherium
addpath(genpath(pwd),'-begin')
spherium_gui

%End of code
function[]=get_nodesets()

% integer mask for set type:
% 1 --> input
% 2 --> measurements
% 3 --> control

% integer mask for operator:
% 1 --> g() <
% 2 --> g() <=
% 3 --> g() ==
% 4 --> g() >=
% 5 --> g() >
% 6 --> < g() >
% 7 --> <= g() >
% 8 --> < g() >=
% 9 --> <= g() >=
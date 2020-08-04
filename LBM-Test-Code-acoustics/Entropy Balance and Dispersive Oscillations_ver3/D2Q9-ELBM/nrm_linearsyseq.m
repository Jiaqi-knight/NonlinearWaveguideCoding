%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% * Author: Ganindu Nanayakkara *                             %%
%  * Asian Institute of Technology *                           %%
% ganindu@gmail.com / facebook.com/ganindu                     %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $$ Free For All $$                                          %%
% *Newton Raphson Method For Linear systems of Equations, demo %%
%                                                              %%
% @paramters: Input Functions: mx, consts, guess               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is specially to solve linear systems of equations using the newton raphson method

%% input data for [mx][x] = [consts] type data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  example: 3 * x1 + 4 * x2 + 5 * x3 = 49 ------ (1)                     %%
%          4 * x1 + 7 * x2 + 3 * x3 = 51 ------ (2)                      %%
%          5 * x1 + 6 * x2 + 5 * x3 = 61 ------ (3)                      %%
% eg:                                                                    %%
%     mx = [ 3 4 5 ;     consts = [49;51;61];    guess  = [5 ;5; 5];     %%
%            4 7 3 ;                                                     %%
%           5 6 5 ];                                                     %%
%                                                                        %%
% user has to apply the values for mx, const and guess                   %%
% the system assigns the x1, x2, ......xn variables and sloves           %%
% the system of equations the input mx, consts and guess matrices        %%
% can be easily fed in using simple file i/o commands.                   %%
% Maximum error and maximum number of iterations can be set also, the    %%
% default values are maxerr = 0.001 and maxit = 10;                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear;clc;
format short g;

mx = [ 3 4 5 ;
       4 7 3 ;
       5 6 5 ]; %matrix for Q1

consts = [49;51;61];
guess  = [5 ;5; 5];

%% Process setup

x = sym('x',[length(mx) 1]);
m = mx * x;
m = m - consts;
J = jacobian(m,x);

iterations  = 0;
maxerr = 0.001;      % Maximum error setting.
maxit = 10;          % Maximum number of iterations setting.
herrx = inf;
current_x = subs(x, x,guess);
%% Iteration
while herrx > maxerr  && iterations < maxit
      
    iterations =  iterations + 1;
    
    c = subs(m,x,current_x);
    xelta= (J^-1)* c;   % can use the "\" also
    xref = current_x;
    partx = eval(-xelta);
    current_x = current_x + partx;
    herrx = max(abs((current_x - xref) ./ current_x)* 100); %highest error precentage
    
    err_log(:,iterations) = herrx;
    result_log(:,iterations) = current_x;
end
%% Results
%% Results matrix will have the form: [iteration, x1, x2, x3,.....,xn, error ] 

master_log = [1:iterations; result_log; err_log]';
disp(master_log);
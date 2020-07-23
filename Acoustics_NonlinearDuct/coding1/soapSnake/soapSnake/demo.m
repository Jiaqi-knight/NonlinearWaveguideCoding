close all
clear all
clc

%% Example 1)
% Euler's formula exp(1i*x)
% Introduction: the most known complex function
   x = linspace(-pi,pi,201); 
   z = exp(1i*x);
   figureFULL; soapSnake(x,z);
   

%% Example 2) 
% The root function ?(x) and its branches 
% To demonstrate how two branches can be displayed (separated with 1/0)
    x = [linspace(-3,0,100), 1/0, linspace(0,+3,100)];
    z =  sqrt(x);
    xInput = [x, 1/0, x]; 
    zInput = [z, 1/0, -z]; 
    figureFULL; soapSnake(xInput,zInput);
    
    
%% Example 3) 
% The (complex) log-function
% To demonstrate that the script automatically detects a isinf-number
% and separates the function in two parts: 
      x = linspace(-6,6,301); % the 51 index will give a problem
      z = log(x); 
      figureFULL; soapSnake(x,z);
      
      
%% Example 4) 
% Here the function automatically divides the function in two branches
      x = linspace(-5,1.6,301);
      z = x.^x; 
      figureFULL; soapSnake(x,z);
      
     
%% Example 5) 
% Here the function automatically divides the function in two branches
    x = linspace(-pi,pi,201); 
    z = exp(1i*2*pi*x).*exp(-1/2*x.^2);
    figureFULL; soapSnake(x,z);
   

    
      
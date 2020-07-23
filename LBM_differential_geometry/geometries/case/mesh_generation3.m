clear all
close all
clc
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
% subfunction_path4=genpath('C:\Users\wjq\Desktop\workspace\geometry-master\geometry')
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
% addpath(subfunction_path4);

formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 50;

x1 = 0;
x2 = 0;
x3_4 = 12;

c = 5;

airfoil = NACAxxxx(Nx,x1,x2,x3_4,c);

f = figure;
plot(airfoil(:,1),airfoil(:,2),'.')
hold on
grid on
xlabel('x')
ylabel('y')
title(strcat('NACA',num2str(x1),num2str(x2),num2str(x3_4),' profile'))
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 50;
Ny = 60;

L = 10;

croot = 5;
ctip  = 2;

x1func = @(y) 0;
x2func = @(y) 0;
x3_4func = @(y) 12;

cfunc = @(y,ctip,croot,L) y*(ctip-croot)/L + croot;

meanlinefunc = @(y) [0.5*y y 0];

[wing,halfwing,upperhalfwing,lowerhalfwing] = wing(Nx,Ny,L,x1func,x2func,x3_4func,cfunc,ctip,croot,meanlinefunc);

f = figure;
plot3(wing(:,1),wing(:,2),wing(:,3),'.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('3D wing')
axis equal

f = figure;
plot3(halfwing(:,1),halfwing(:,2),halfwing(:,3),'.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('3D wing')
axis equal

f = figure;
plot3(upperhalfwing(:,1),upperhalfwing(:,2),upperhalfwing(:,3),'.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('3D wing')
axis equal

f = figure;
plot3(lowerhalfwing(:,1),lowerhalfwing(:,2),lowerhalfwing(:,3),'.')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('3D wing')
axis equal
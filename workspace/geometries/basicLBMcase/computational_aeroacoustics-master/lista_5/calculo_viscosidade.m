clc; clear; close all;
% Calculo da viscosidade para achar a discretizacao x
viscosidade_ar = 1.5e-5; % m^2/s
velocidade_som_fisica = 343; % m/s
omega = 1.93; % Relaxation frequency
tau = 1/omega;
A = viscosidade_ar*sqrt(3);
B = velocidade_som_fisica*(tau - 1/2);
Dx = A/B;
fprintf('Tamanho do delta X: %f milimetros. \n', 1000*Dx);
fprintf('Dimensoes fisicas do anteparo: %f milimetros por %f milimetros.\n', 1000*Dx*40, 1000*Dx*40);
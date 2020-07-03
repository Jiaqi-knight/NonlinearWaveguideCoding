% Script para calcular as pressoes no campo distante
%close all; clear('all'); clc;
%load('dados.mat');
tempo_total = 1000;
frequencia = (0.0001077+0.02*tempo_total)/tempo_total;
% y e x
delta_x = 0.004177e-3; % metros
%posicao_observador = [15/delta_x 0]; % metros convertidos para lattice
cs = 1/sqrt(3);
cs2 = cs^2;
Ma = [0.03 0.07 0.1];
U1 = 0;%Ma(3)*cs;
U2 = 0;
ponto_1_superficie = [(231 - 1) (231 - 1)];
ponto_2_superficie = [(271 + 1) (271 + 1)];

theta = (0:1:180)*pi/180;
%theta = 0:0.02:2*pi;
%Angle = 0;
raio = 15/delta_x;
pressoes = calcular_campo_distante(ponto_1_superficie, ... 
ponto_2_superficie, U1, U2, frequencia, theta, raio);
figure; polar(theta, pressoes);



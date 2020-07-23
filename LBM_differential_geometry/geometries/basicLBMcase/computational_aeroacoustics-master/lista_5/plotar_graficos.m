clear('all'); clc; close all;
% plotar graficos dos resultados
mach_007_1 = load('dados_results_M_0.07_1.mat');
mach_007_50 = load('dados_results_M_0.07_50.mat');
mach_007_100 = load('dados_results_M_0.07_100.mat');
mach_01_1 = load('dados_results_M_0.1_1.mat');
mach_01_50 = load('dados_results_M_0.1_50.mat');
mach_01_100 = load('dados_results_M_0.1_100.mat');
theta = mach_007_1.data(1, :);

% plotando somente os dipolos, variando numero de celulas
figure(1); 
a = abs(mach_007_50.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'red');
hold on;
a = abs(mach_007_1.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'blue'); hold on;
a = abs(mach_007_100.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'black');
title('Dipolo - Mach 0.07');
legend('50 celulas', '1 celula', '100 celulas');

a = abs(mach_01_1.data(3, :));
a = (a.^2)/mean(a)^2;
figure(2); polar(theta, a, 'blue'); 
hold on;
a = abs(mach_01_50.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'red');
hold on;
a = abs(mach_01_100.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'black');
title('Dipolo - Mach 0.1');
legend('1 celula', '50 celulas', '100 celulas');

% plotando somente os dipolos, variando o Mach
figure(3); 
a = abs(mach_01_1.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'red');
hold on;
a = abs(mach_007_1.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'blue'); 
title('Dipolo - 1 celula');
legend('Mach 0.1', 'Mach 0.07');

figure(4); 
a = abs(mach_01_50.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'red');
hold on;
a = abs(mach_007_50.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'blue'); 
title('Dipolo - 50 celulas');
legend('Mach 0.1', 'Mach 0.07');

figure(5);
a = abs(mach_01_100.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'red');
hold on;
a = abs(mach_007_100.data(3, :));
a = (a.^2)/mean(a)^2;
polar(theta, a, 'blue'); 
title('Dipolo - 100 celulas');
legend('Mach 0.1', 'Mach 0.07');

% plotando grafico das intensidades frequencias pela posicao
figure(6);
load('F_i_peak.mat');
F_i_peak_1 = F_i_peak(1:42);
F_i_peak_2 = F_i_peak(1*42:2*42-1);
F_i_peak_3 = F_i_peak(2*42:3*42-1);
F_i_peak_4 = F_i_peak(3*42:4*42-1);
plot([1:length(F_i_peak_1)], F_i_peak_1, 'blue');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_2, 'red');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_3, 'black');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_4, 'cyan');
title('Pressao do Dipolo ao Longo da Superficie - 1 celula');
xlabel('Numero de Celulas de Lattice');
ylabel('Pressao');
legend('Lado 1', 'Lado 2', 'Lado 3', 'Lado 4');

figure(7);
load('Q_peak.mat');
Q_peak_1 = Q_peak(1:42);
Q_peak_2 = Q_peak(1*42:2*42-1);
Q_peak_3 = Q_peak(2*42:3*42-1);
Q_peak_4 = Q_peak(3*42:4*42-1);
plot([1:length(Q_peak_1)], Q_peak_1, 'blue');
hold on;
plot([1:length(Q_peak_1)], Q_peak_2, 'red');
hold on;
plot([1:length(Q_peak_1)], Q_peak_3, 'black');
hold on;
plot([1:length(Q_peak_1)], Q_peak_4, 'cyan');
title('Pressao do Monopolo ao Longo da Superficie - 1 celula');
xlabel('Numero de Celulas de Lattice');
ylabel('Pressao');
legend('Lado 1', 'Lado 2', 'Lado 3', 'Lado 4');

% plotando grafico das intensidades frequencias pela posicao
figure(8);
load('Fi_peak_M_0.1_100.mat');
F_i_peak_1 = Fi_peak_100(1:240);
F_i_peak_2 = Fi_peak_100(1*240:2*240-1);
F_i_peak_3 = Fi_peak_100(2*240:3*240-1);
F_i_peak_4 = Fi_peak_100(3*240:4*240-1);
plot([1:length(F_i_peak_1)], F_i_peak_1, 'blue');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_2, 'red');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_3, 'black');
hold on;
plot([1:length(F_i_peak_1)], F_i_peak_4, 'cyan');
title('Pressao do Dipolo ao Longo da Superficie - 100 celulas');
xlabel('Numero de Celulas de Lattice');
ylabel('Pressao');
legend('Lado 1', 'Lado 2', 'Lado 3', 'Lado 4');

figure(9);
load('Q_peak_M_0.1_100.mat');
Q_peak_1 = Q_peak_100(1:240);
Q_peak_2 = Q_peak_100(1*240:2*240-1);
Q_peak_3 = Q_peak_100(2*240:3*240-1);
Q_peak_4 = Q_peak_100(3*240:4*240-1);
plot([1:length(Q_peak_1)], Q_peak_1, 'blue');
hold on;
plot([1:length(Q_peak_1)], Q_peak_2, 'red');
hold on;
plot([1:length(Q_peak_1)], Q_peak_3, 'black');
hold on;
plot([1:length(Q_peak_1)], Q_peak_4, 'cyan');
title('Pressao do Monopolo ao Longo da Superficie - 100 celulas');
xlabel('Numero de Celulas de Lattice');
ylabel('Pressao');
legend('Lado 1', 'Lado 2', 'Lado 3', 'Lado 4');

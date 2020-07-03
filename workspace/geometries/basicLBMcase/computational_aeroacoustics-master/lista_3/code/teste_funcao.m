% script de teste da funcao 'funcao'
close all;
clear('all');

i = 1;
vecx = [1:40] + 0.3;
vecy = 20*sin([1:40] + 0.3) + 21;
%vecy = ([1:40] + 0.3);
Nr = 50;
Mc = 50;
[q, vecloc] = funcao(i, vecx, vecy, Nr, Mc);

figure;
h1 = axes;
imagesc(vecloc);
hold on;
plot(vecx, vecy, 'black');
grid on;
set(h1, 'Ydir', 'normal');
legend('Parede Nao-Alinhada');

figure;
h1 = axes;
imagesc(q);
hold on;
plot(vecx, vecy, 'black');
grid on;
set(h1, 'Ydir', 'normal');
legend('Parede Nao-Alinhada');

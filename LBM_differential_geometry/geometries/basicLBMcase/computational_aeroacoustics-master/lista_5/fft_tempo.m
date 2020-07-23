clc; clear('all'); close all;
cs = 1/sqrt(3);
cs2 = cs^2;

% fazer o load aqui do tempo total
load('dados/1.mat');
tempo_total = total_time - 15000;
pressoes(1:tempo_total) = 0;
for ta = 1:tempo_total
    if  mod(ta, 1000) == 0
        disp('Progresso para os Fi: ');
        disp((ta*100)/tempo_total);
    end
    nome_arquivo = ['dados/' num2str(ta) '.mat'];
    load(nome_arquivo);
    % Abrindo a matriz rho para pegar a densidade naquele ponto no time-step
    linha_3 = rho_linha(2*43 + 1 : 3*43);

    % Calculando a pressao naquele ponto num determinado time-step
    pressao = linha_3((43 - 1)/2 + 1)*cs2

    % Inserindo o novo ponto no conjunto de pontos no tempo
    pressoes(tempo_total) = pressao;
end

% Realizando a FFT
% Plotando o valor dela
frequencias = linspace(0, 1, tempo_total);
figure; plot(frequencias, abs(fft(pressoes)));
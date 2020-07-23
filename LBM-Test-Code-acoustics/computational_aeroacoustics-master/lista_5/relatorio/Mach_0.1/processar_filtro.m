clear('all');
close all;
clc;

% capturar campo de pressao
load rho;
campo_acustico = (rho - 1)/3;

% transformada de fourier para cada linha
fft_campo_acustico = fft(campo_acustico, length(campo_acustico), 2);

% filtro a transformada de fourier usando a hanning dado o numero de onda k
strouhal_grafico = 0.08164;
mach = 0.1;
velocidade_som = 1/sqrt(3);
dimensao_cubo = 40;
freq_pico = (strouhal_grafico*mach*velocidade_som)/dimensao_cubo;
k = 2*pi*freq_pico/(velocidade_som);
frequencias = linspace(0, 1, length(fft_campo_acustico(350, 30:502-30)));
%strouhals = frequencias*dimensao_cubo/mach*velocidade_som;
k_s = 2*pi*frequencias/(velocidade_som);
figure; plot(k_s, abs(fft_campo_acustico(400, 30:502-30)), '*');
%hold on;
%plot(imag(fft_campo_acustico(1, :)), 'r');



% transformada inversa linha por linha 
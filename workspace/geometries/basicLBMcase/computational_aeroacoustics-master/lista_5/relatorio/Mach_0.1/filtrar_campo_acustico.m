clear('all'); clc;% close all;

% Filtrando campo acustico
%% Calculando o campo acustico
load rho.mat
[linhas colunas] = size(rho);
velocidade_som_lattice = 1/sqrt(3);
campo_acustico = (rho - 1)*(velocidade_som_lattice^2);
campo_acustico = campo_acustico/max(max(campo_acustico));

%% Calculando o filtro para convoluir
numero_strouhal = 0.08164;
frequencia = numero_strouhal*(0.1*(1/sqrt(3)))/40;
numero_onda = 2*pi*frequencia/(1/sqrt(3));
times = 0 : colunas - 1;
filtro = chirp(times, ... 
numero_onda, times(end)*2, numero_onda);

%% Calculando campo acustico filtrado total
campo_acustico_filtrado = conv2(filtro, filtro, campo_acustico);
campo_acustico_filtrado = campo_acustico_filtrado/max(max(campo_acustico_filtrado));
figure;
subplot(1,2,1);
imagesc(campo_acustico(1:linhas, 1:colunas), [-.15 .15]);
axis equal;
subplot(1,2,2);
imagesc(campo_acustico_filtrado(1:linhas, 1:colunas), [-.5 .5]);
axis equal;
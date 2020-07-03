% struct pressure
load struct_pressure;

Angle = (0:3:180)*pi/180;
diretividade(1:length(Angle)) = 0;
diretividades = [];
for angulo = 1:length(Angle)
	vetor_pressao = struct_pressure{angulo};
	fft_pressao = abs(fft(vetor_pressao(14500:end)));
	diretividade(angulo) = fft_pressao(8);
end

figure; polar(Angle, diretividade);
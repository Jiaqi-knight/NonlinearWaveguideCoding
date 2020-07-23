%% build_anechoic_condition: function description
function [sigma_mat9 Ft] = build_anechoic_condition(number_lines_lattice, ... 
number_columns_lattice, distance, growth_delta)

	lattice_sound_speed = 1/sqrt(3);
	lattice_sound_speed_pow_2 = lattice_sound_speed^2;
	w1=4/9;     % centro    %pesos de relaxaçao devido ao D2Q9 (pg.20)
	w2=1/9;     % ortogonais        
	w3=1/36;    % diagonais
	coef1=  1/(2*lattice_sound_speed_pow_2^2); %para uso na relaxacao
	coef2= -1/(2*lattice_sound_speed_pow_2);
	% funcoes distribuiçao (Eq. 1.46)
	Ft = zeros(number_lines_lattice,number_columns_lattice,9);
	%funcoes target
	Ux_t=0;
	Uy_t=0;
	U_t=Ux_t^2+Uy_t^2;
	densi_t = 1;
	Ft(:,:,9)= w1*densi_t.*(1+coef2*U_t);
	Ft(:,:,1)= w2*densi_t.*(1 +Ux_t/lattice_sound_speed_pow_2 +coef1*(Ux_t.^2 )+coef2*U_t);
	Ft(:,:,2)= w2*densi_t.*(1 +Uy_t/lattice_sound_speed_pow_2 +coef1*(Uy_t.^2 ) +coef2*U_t);
	Ft(:,:,3)= w2*densi_t.*(1 -Ux_t/lattice_sound_speed_pow_2 +coef1*(Ux_t.^2 )+coef2*U_t);
	Ft(:,:,4)= w2*densi_t.*(1 -Uy_t/lattice_sound_speed_pow_2 +coef1*(Uy_t.^2) +coef2*U_t);
	Ft(:,:,5)= w3*densi_t.*(1 +(+Ux_t+Uy_t)/lattice_sound_speed_pow_2 +coef1*((+Ux_t+Uy_t).^2) +coef2*U_t);
	Ft(:,:,6)= w3*densi_t.*(1 +(-Ux_t+Uy_t)/lattice_sound_speed_pow_2 +coef1*((-Ux_t+Uy_t).^2) +coef2*U_t);
	Ft(:,:,7)= w3*densi_t.*(1 +(-Ux_t-Uy_t)/lattice_sound_speed_pow_2 +coef1*((-Ux_t-Uy_t).^2) +coef2*U_t);
	Ft(:,:,8)= w3*densi_t.*(1 +(+Ux_t-Uy_t)/lattice_sound_speed_pow_2 +coef1*((+Ux_t-Uy_t).^2) +coef2*U_t);
	% 
	D_t = distance;  % em número de celulas
	sigma_t = 0.3;
	delta_t = 0:D_t;
	sigma = sigma_t*(delta_t/D_t).^2;
	sigma_mat = [];
	for i = 1:number_lines_lattice  % ver se tem jeito melhor de concatenar as matrizes
	    sigma_mat = cat(1,sigma,sigma_mat);
	end
	if growth_delta == 1
		sigma_mat = [zeros(number_lines_lattice,number_columns_lattice-D_t-1) sigma_mat];
	elseif growth_delta == -1
		sigma_mat = fliplr([zeros(number_lines_lattice,number_columns_lattice-D_t-1) sigma_mat]);
	else
		sigma_mat = [zeros(number_lines_lattice,number_columns_lattice-D_t-1) sigma_mat];
	end
	

	sigma_mat9 = [];
	for i = 1:9
	sigma_mat9 = cat(3,sigma_mat,sigma_mat9);
	end
	% Condicao anecoica no final

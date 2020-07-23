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
	Ma=[0.03 0.07 0.1];
	Ux_t=Ma(3)*lattice_sound_speed;
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
	delta_t = 0:D_t - 1;
	sigma = sigma_t*(delta_t/D_t).^2;
	sigma = sigma';
	sigma_mat = [];

	% x e y nesse caso
	sigma_mat9 = zeros(number_lines_lattice, number_columns_lattice, 9);
	matrix_sigma_leaf = zeros(number_lines_lattice, number_columns_lattice);
	% condicao anecoica para a direita
	if growth_delta == 1
		% construindo a matriz de nove folhas sigma
		size_x = number_lines_lattice - D_t + 1 : number_lines_lattice;
		size_y = 1:number_columns_lattice;
		matrix_sigma_leaf(size_x, size_y) = 1;
		
		for	number_column = 1: number_columns_lattice
			matrix_sigma_leaf(size_x, number_column) = ...
			matrix_sigma_leaf(size_x, number_column).*sigma;
		end
	% condicao anecoica para a esquerda
	elseif growth_delta == -1
		% construindo a matriz de nove folhas sigma
		size_x = 1 : D_t;
		size_y = 1 : number_columns_lattice;
		matrix_sigma_leaf(size_x, size_y) = 1;
		for	number_column = 1: number_columns_lattice
			matrix_sigma_leaf(size_x, number_column) = ...
			matrix_sigma_leaf(size_x, number_column).*flip(sigma);
		end
	% condicao anecoica para cima
	elseif growth_delta == 0.5
		size_x = 1 : number_lines_lattice;
		size_y = number_columns_lattice - D_t + 1 : number_columns_lattice;
		matrix_sigma_leaf(size_x, size_y) = 1;
		for number_line =  1 : number_lines_lattice
			matrix_sigma_leaf(number_line, size_y) = ...
			matrix_sigma_leaf(number_line, size_y).*sigma';
		end
	% condicao anecoica para baixo
	elseif growth_delta == -0.5
		size_x = 1 : number_lines_lattice;
		size_y = 1 : D_t;
		matrix_sigma_leaf(size_x, size_y) = 1;
		for number_line =  1 : number_lines_lattice
			matrix_sigma_leaf(number_line, size_y) = ...
			matrix_sigma_leaf(number_line, size_y).*flip(sigma)';
		end
	end
	
	for number_leaf = 1:9
		sigma_mat9(:,:,number_leaf) = matrix_sigma_leaf;
	end

	% Agora eh y por x
	sigma_mat9_aux = zeros(number_columns_lattice, number_lines_lattice, 9);
	Ft_aux = zeros(number_columns_lattice, number_lines_lattice, 9);
	for direction = 1 : 9
		sigma_mat9_aux(:, :, direction) = sigma_mat9(:, :, direction)';
		Ft_aux(:, :, direction) = Ft(:, :, direction)';
	end 
	sigma_mat9 = sigma_mat9_aux;
	Ft = Ft_aux;

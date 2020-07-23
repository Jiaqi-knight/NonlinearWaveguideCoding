%% build_anechoic_condition: function description
function [sigma_mat9 Ft] = build_anechoic_condition_axis(number_lines_lattice, ... 
number_columns_lattice, distance, growth_delta, e, e_alpha, w_alpha)

	% funcoes distribuicao (Eq. 1.46)
	Ft = zeros(number_lines_lattice,number_columns_lattice,9);
	%funcoes target
	ux=0;
	uy=0;
	densi_t = 1;
	
    for link = 1:9
        c1 = 3/(e^2);
        C1 = e_alpha(link,2)*ux + e_alpha(link,1)*uy;
        c2 = 9/(2*e^4);
        C2 = (e_alpha(link,2)^2)*(ux.^2) + 2*e_alpha(link,1)*e_alpha(link,2)*uy.*ux ...
        + (e_alpha(link,1)^2)*(uy.^2);
        c3 = 3/(2*e^2);
        C3 = ux.^2 + uy.^2;
        Ft(:,:,link)= w_alpha(link)*densi_t .*(1 + c1*C1 + c2*C2 - c3*C3);
        %mean(mean(feq(:,:,link)))
        %link
    end

	% 
	D_t = distance;  % em numero de celulas
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

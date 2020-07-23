%% build_source_anechoic: direction = 1 (virado para cima) 
%						  direction = -1 (virado para baixo)
%						  direction = 0.5 (virado para direita)
%						  direction = -0.5 (virado para esquerda)
function [sigma_source Ft_source] = build_source_anechoic(Nr, Mc, ...
density, Ux_t, Uy_t, point_y, point_x, distance_x, distance_y, direction)
	
	lattice_sound_speed = 1/sqrt(3);
	lattice_sound_speed_pow_2 = lattice_sound_speed^2;
	cs2 = lattice_sound_speed_pow_2;
	cs = lattice_sound_speed;
	w1=4/9;     % centro    %pesos de relaxaçao devido ao D2Q9 (pg.20)
	w2=1/9;     % ortogonais        
	w3=1/36;    % diagonais
	coef1=  1/(2*cs2^2); %para uso na relaxacao
	coef2= -1/(2*cs2);
	% funcoes distribuiçao (Eq. 1.46)
	Ft_source = zeros(Nr, Mc, 9);
	%funcoes target
	Ux_t = Ux_t;
	Uy_t = Uy_t;
	U_t=Ux_t^2+Uy_t^2;
	densi_t = density;
	Ft_source(:,:,9)= w1*densi_t.*(1+coef2*U_t);
	Ft_source(:,:,1)= w2*densi_t.*(1 +Ux_t/cs2 +coef1*(Ux_t.^2 )+coef2*U_t);
	Ft_source(:,:,2)= w2*densi_t.*(1 +Uy_t/cs2 +coef1*(Uy_t.^2 ) +coef2*U_t);
	Ft_source(:,:,3)= w2*densi_t.*(1 -Ux_t/cs2 +coef1*(Ux_t.^2 )+coef2*U_t);
	Ft_source(:,:,4)= w2*densi_t.*(1 -Uy_t/cs2 +coef1*(Uy_t.^2) +coef2*U_t);
	Ft_source(:,:,5)= w3*densi_t.*(1 +(+Ux_t+Uy_t)/cs2 +coef1*((+Ux_t+Uy_t).^2) +coef2*U_t);
	Ft_source(:,:,6)= w3*densi_t.*(1 +(-Ux_t+Uy_t)/cs2 +coef1*((-Ux_t+Uy_t).^2) +coef2*U_t);
	Ft_source(:,:,7)= w3*densi_t.*(1 +(-Ux_t-Uy_t)/cs2 +coef1*((-Ux_t-Uy_t).^2) +coef2*U_t);
	Ft_source(:,:,8)= w3*densi_t.*(1 +(+Ux_t-Uy_t)/cs2 +coef1*((+Ux_t-Uy_t).^2) +coef2*U_t);

	% posicionando a fonte
	sigma_source = zeros(Nr, Mc, 9);
	% assintota para a direita
	if direction == 0.5
		sigma_t = 0.3;
		delta_t = 0:distance_x;
		sigma = sigma_t*(delta_t/distance_x).^2;
		sigma_source_leaf = zeros(Nr, Mc);
		size_x = point_x : point_x + distance_x;
		size_y = point_y : point_y + distance_y;
		sigma_source_leaf(size_y, size_x) = 1;
		for point_y = size_y(1) : size_y(end)
			sigma_source_leaf(point_y, size_x) = ... 
			sigma_source_leaf(point_y, size_x).*flip(sigma);
		end
		for direction = 1:9
			sigma_source(:, :, direction) = sigma_source_leaf(:, :);
		end
	end
	
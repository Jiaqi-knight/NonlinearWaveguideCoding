function [q, vecloc] = funcao(i, vecx, vecy, Nr, Mc)

	% verify if the points is inside of lattice
	is_not_points_inside_lattice = ... 
	min(vecx) < 1 || max(vecx) > Mc ...
	|| min(vecy) < 1 || max(vecy) > Nr;
	if is_not_points_inside_lattice
		i = -1;
		message = 'The points x and y are outside of lattice.';
		disp(message);
	end

	% direction 1 of lattice cell 
	if i == 1
		vecloc = zeros(Nr, Mc);
		q = zeros(Nr, Mc);
		% square that have the entire points
		square_x = [floor(min(vecx)) ceil(max(vecx))];
		square_y = [floor(min(vecy)) ceil(max(vecy))];
		%vecloc(square_x(1):square_x(2), ...
		%square_y(1):square_y(2)) = 1;
		
		% getting points on the left of curve
		distances_points = [];
		for point_y = square_y(1):square_y(2)
			% looking for a x point 
			% (Verify whats p1 and p2 is nearest from height)
			distances_y = abs(vecy - point_y);
			slot_min = find(distances_y == min(distances_y));
			p2_y = vecy(slot_min);
			p2_x = vecx(slot_min);
			distances_y(slot_min) = 10e10;
			slot_min = find(distances_y == min(distances_y));
			p1_y = vecy(slot_min);
			p1_x = vecx(slot_min);

			% Agora tenho p1 e p2, tenho que agora fazer a
			% interpolacao linear para achar o p3 
			% equacao da reta: y = ax + b
			a = (p2_y - p1_y)/(p2_x - p1_x);
			b = p1_y - a*p1_x;
			p3_x = (point_y - b)/a;
			
			% e calcular a distancia em x do ponto que eu to para o x de p3
			distances_points = [square_x(1):square_x(2)] - p3_x;
			distances_points = abs(distances_points);
			point_x = find(distances_points == min(distances_points));
			vecloc(point_y, point_x) = 1;
			q(point_y, point_x) = min(distances_points);
		end


	else
		q = 0; vecloc = 0;
	end
	
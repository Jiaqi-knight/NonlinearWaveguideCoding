%% build_source_anechoic: direction = 1 (virado para cima) 
%             direction = -1 (virado para baixo)
%             direction = 0.5 (virado para direita)
%             direction = -0.5 (virado para esquerda)
function [sigma_source Ft] = build_source_anechoic_axis(Nr, Mc, ...
density_source, Ux_t, Uy_t, point_y, point_x, distance_x, distance_y, direction, ...
e, e_alpha, w_alpha)
  
  % funcoes distribuicao (Eq. 1.46)
  Ft = zeros(Nr,Mc,9);
  %funcoes target
  ux=Ux_t;
  uy=Uy_t;
  densi_t = density_source;
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
  
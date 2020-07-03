%% calcular_campo_distante: para cada time-step eh um arquivo na pasta 'dados'
function pressoes = calcular_campo_distante(ponto_1_superficie, ponto_2_superficie, ...
U1, U2, frequencia, Angle, raio)
cs = 1/sqrt(3);
cs2 = cs^2;

% fazer o load aqui do tempo total
load('dados/1.mat');
tempo_total = total_time - 42000;

% Efetuando os calculos de F_i
F_i = {};
for ta = 1:tempo_total
    if  mod(ta, 1000) == 0
        disp('Progresso para os Fi: ');
        disp((ta*100)/tempo_total);
    end
    
    nome_arquivo = ['dados/' num2str(ta) '.mat'];
    load(nome_arquivo);

    % Montando os planos
    % invertendo rho linha
    tamanho_cada_linha = ponto_2_superficie(1) - ponto_1_superficie(1) + 1;
    rho(1:tamanho_cada_linha, 1:tamanho_cada_linha) = 0;
    y1 = 1;
    x1 = 1;
    y2 = tamanho_cada_linha;
    x2 = tamanho_cada_linha;
    rho(y1:y1, x1:x2) = rho_linha(1:tamanho_cada_linha);
    rho(y1:y2, x2:x2) = rho_linha(tamanho_cada_linha + 1:2*tamanho_cada_linha)';
    rho(y2:y2, x2:-1:x1) = rho_linha(2*tamanho_cada_linha + 1:3*tamanho_cada_linha);
    rho(y2:-1:y1,x1:x1) = rho_linha(3*tamanho_cada_linha + 1:4*tamanho_cada_linha)';
    % invertendo ux
    ux(1:tamanho_cada_linha, 1:tamanho_cada_linha) = 0;
    ux(y1:y1, x1:x2) = ux_linha(1:tamanho_cada_linha);
    ux(y1:y2, x2:x2) = ux_linha(tamanho_cada_linha + 1:2*tamanho_cada_linha)';
    ux(y2:y2, x2:-1:x1) = ux_linha(2*tamanho_cada_linha + 1:3*tamanho_cada_linha);
    ux(y2:-1:y1,x1:x1) = ux_linha(3*tamanho_cada_linha + 1:4*tamanho_cada_linha)';
    % invertendo uy
    uy(1:tamanho_cada_linha, 1:tamanho_cada_linha) = 0;
    uy(y1:y1, x1:x2) = uy_linha(1:tamanho_cada_linha);
    uy(y1:y2, x2:x2) = uy_linha(tamanho_cada_linha + 1:2*tamanho_cada_linha)';
    uy(y2:y2, x2:-1:x1) = uy_linha(2*tamanho_cada_linha + 1:3*tamanho_cada_linha);
    uy(y2:-1:y1,x1:x1) = uy_linha(3*tamanho_cada_linha + 1:4*tamanho_cada_linha)';

    % definindo velocidades
    u1 = ux;
    u2 = uy;
    rho_l = 1;

    Fi = ux;
    Fi(:) = 0;
    % Para baixo ni(y x) = [-1 0] (posicao_y == y1 && posicao_x <= x2)
    % Para direita ni(y x) = [0 1] (posicao_y <= y2 && posicao_x == x2)
    % Para cima ni(y x) = [1 0] (posicao_y == y2 && posicao_x <= x2)
    % Para esquerda ni(y x) = [0 -1] (posicao_y <= y2 && posicao_x == x1)
    for posicao_y = 1:tamanho_cada_linha
        for posicao_x = 1:tamanho_cada_linha
            % baixo
            if posicao_y == y1 && posicao_x <= x2 
                n1 = 0;
                n2 = -1;

                %
                p = (rho(posicao_y, posicao_x) - 1)*cs2;
                rho_u1_2U1_u1 = rho(posicao_y, posicao_x).*u1(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1);
                rho0_U1_2 = rho_l*U1^2;
                F1_n1 = p + rho_u1_2U1_u1 + rho0_U1_2;
                %
                rho_u1_2U1_u2 = rho(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1).*u2(posicao_y, posicao_x);
                rho0_U1_U2 = rho_l*U1*U2;
                F1_n2 = rho_u1_2U1_u2 + rho0_U1_U2;
                %
                F1 = F1_n1*n1 + F1_n2*n2;
                %
                rho_u2_2U2_u1 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u1(posicao_y, posicao_x);
                F2_n1 = rho_u2_2U2_u1 + rho0_U1_U2;
                %
                rho_u2_2U2_u2 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u2(posicao_y, posicao_x);
                rho0_U2_2 = rho_l*U2^2;
                F2_n2 = p + rho_u2_2U2_u2 + rho0_U1_U2 + rho0_U2_2; 
                %
                F2 = F2_n1*n1 + F2_n2*n2;
                %
                Fi(posicao_y, posicao_x) = sqrt(F1^2 + F2^2);

            % direita
            elseif posicao_y <= y2 && posicao_x == x2
                n1 = 1;
                n2 = 0;

                %
                p = (rho(posicao_y, posicao_x) - 1)*cs2;
                rho_u1_2U1_u1 = rho(posicao_y, posicao_x).*u1(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1);
                rho0_U1_2 = rho_l*U1^2;
                F1_n1 = p + rho_u1_2U1_u1 + rho0_U1_2;
                %
                rho_u1_2U1_u2 = rho(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1).*u2(posicao_y, posicao_x);
                rho0_U1_U2 = rho_l*U1*U2;
                F1_n2 = rho_u1_2U1_u2 + rho0_U1_U2;
                %
                F1 = F1_n1*n1 + F1_n2*n2;
                %
                rho_u2_2U2_u1 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u1(posicao_y, posicao_x);
                F2_n1 = rho_u2_2U2_u1 + rho0_U1_U2;
                %
                rho_u2_2U2_u2 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u2(posicao_y, posicao_x);
                rho0_U2_2 = rho_l*U2^2;
                F2_n2 = p + rho_u2_2U2_u2 + rho0_U1_U2 + rho0_U2_2; 
                %
                F2 = F2_n1*n1 + F2_n2*n2;
                %
                Fi(posicao_y, posicao_x) = sqrt(F1^2 + F2^2);

            % cima
            elseif posicao_y == y2 && posicao_x <= x2
                n1 = 0;
                n2 = 1;

                %
                p = (rho(posicao_y, posicao_x) - 1)*cs2;
                rho_u1_2U1_u1 = rho(posicao_y, posicao_x).*u1(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1);
                rho0_U1_2 = rho_l*U1^2;
                F1_n1 = p + rho_u1_2U1_u1 + rho0_U1_2;
                %
                rho_u1_2U1_u2 = rho(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1).*u2(posicao_y, posicao_x);
                rho0_U1_U2 = rho_l*U1*U2;
                F1_n2 = rho_u1_2U1_u2 + rho0_U1_U2;
                %
                F1 = F1_n1*n1 + F1_n2*n2;
                %
                rho_u2_2U2_u1 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u1(posicao_y, posicao_x);
                F2_n1 = rho_u2_2U2_u1 + rho0_U1_U2;
                %
                rho_u2_2U2_u2 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u2(posicao_y, posicao_x);
                rho0_U2_2 = rho_l*U2^2;
                F2_n2 = p + rho_u2_2U2_u2 + rho0_U1_U2 + rho0_U2_2; 
                %
                F2 = F2_n1*n1 + F2_n2*n2;
                %
                Fi(posicao_y, posicao_x) = sqrt(F1^2 + F2^2);

            % esquerda
            elseif posicao_y <= y2 && posicao_x == x1
                n1 = -1;
                n2 = 0;

                %
                p = (rho(posicao_y, posicao_x) - 1)*cs2;
                rho_u1_2U1_u1 = rho(posicao_y, posicao_x).*u1(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1);
                rho0_U1_2 = rho_l*U1^2;
                F1_n1 = p + rho_u1_2U1_u1 + rho0_U1_2;
                %
                rho_u1_2U1_u2 = rho(posicao_y, posicao_x).*(u1(posicao_y, posicao_x) - 2*U1).*u2(posicao_y, posicao_x);
                rho0_U1_U2 = rho_l*U1*U2;
                F1_n2 = rho_u1_2U1_u2 + rho0_U1_U2;
                %
                F1 = F1_n1*n1 + F1_n2*n2;
                %
                rho_u2_2U2_u1 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u1(posicao_y, posicao_x);
                F2_n1 = rho_u2_2U2_u1 + rho0_U1_U2;
                %
                rho_u2_2U2_u2 = rho(posicao_y, posicao_x).*(u2(posicao_y, posicao_x) - 2*U2).*u2(posicao_y, posicao_x);
                rho0_U2_2 = rho_l*U2^2;
                F2_n2 = p + rho_u2_2U2_u2 + rho0_U1_U2 + rho0_U2_2; 
                %
                F2 = F2_n1*n1 + F2_n2*n2;
                %
                Fi(posicao_y, posicao_x) = sqrt(F1^2 + F2^2);
            end

                
        end
    end
    F_i{ta} = Fi;
    %size(F_i);
end

% Como ja tenho cada plano ao longo do tempo, agora eh fazer a FFT
%% Calculando Fi_fft
Fi_fft = {};
for ponto_x_plano = 1:tamanho_cada_linha
    for ponto_y_plano = 1:tamanho_cada_linha
        if ponto_y_plano == y1 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x2 || ...
            ponto_y_plano == y2 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x1

            % Mergulhando nas matrizes
            vetor_auxiliar(1:tempo_total) = 0;
            for incremento_tempo = 1:tempo_total
                matriz_Fi = F_i{incremento_tempo};
                vetor_auxiliar(incremento_tempo) = ... 
                matriz_Fi(ponto_y_plano, ponto_x_plano);
            end
            fft_linha_mergulhada = fft(vetor_auxiliar);

            % Guardando os vetores fft's
            Fi_fft{ponto_y_plano, ponto_x_plano} = fft_linha_mergulhada;

        else

        end
    end
end

% Desacoplando os vetores de Fi_fft para montar as superficies certinhas
fft_Fi = {};
for incremento_tempo = 1:tempo_total
    matriz_auxiliar = F_i{incremento_tempo};
    matriz_auxiliar(:) = 0;
    for ponto_y_plano = 1:tamanho_cada_linha
        for ponto_x_plano = 1:tamanho_cada_linha
            if ponto_y_plano == y1 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x2 || ...
            ponto_y_plano == y2 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x1

            vetor_fft_auxiliar = Fi_fft{ponto_y_plano, ponto_x_plano};
            matriz_auxiliar(ponto_y_plano, ponto_x_plano) = ...
            vetor_fft_auxiliar(incremento_tempo);

            end
        end
    end
    fft_Fi{incremento_tempo} = matriz_auxiliar;
end

% Selecionando o valor da frequencia
slot_frequencia = 5;
Fi_frequencia_escolhida = fft_Fi{slot_frequencia};
frequencia = slot_frequencia/tempo_total;
%frequencia = 0.0001;

% -------------------------------------------
pressoes = Angle;
pressoes(:) = 0;
for ponto = 1 : length(Angle)
   x = round(raio*cos(Angle(ponto)));
   y = round(raio*sin(Angle(ponto)));
   disp('Progresso: ');
   disp(ponto*100/length(Angle));

    % Construindo a funcao Green
    %% Construindo A e B
    A = F_i{1};
    A(:) = 0;
    B = A;
    for epsilon = 1:tamanho_cada_linha
        for neta = 1:tamanho_cada_linha
            %if epsilon == y1 && neta <= x2 || ...
            %epsilon <= y2 && neta == x2 || ...
            %epsilon == y2 && neta <= x2 || ...
            %epsilon <= y2 && neta == x1

                V = 0;
                U = U1;
                %teta = atan(uy(e)/U);
                %x_ = (x - epsilon)*cos(teta) + (y - neta)*sin(teta);
                x_ = (x - epsilon)^2;
                k = 2*pi*frequencia/cs;
                M = 0.1;
                beta = sqrt((1-M^2));
                A(epsilon, neta) = (i/(4*beta))*exp((i*M*k*x_)/beta^2);
                %y_ = -(x - epsilon)*sin(teta) + (y - neta)*cos(teta);
                y_ = (y - neta)^2;
                %z = (k/(beta^2))*sqrt(x_^2 + (beta^2)*(y_^2));
                %H = besselj(0,z) - 1i*bessely(0,z);
                B(epsilon, neta) = besselh(0, 2, (k/(beta^2))*sqrt(x_^2 + (beta^2)*(y_^2)));
                %B(epsilon, neta) = H;
            %end
        end
    end

    funcao_green = A.*B;
    diff_funcao_green = diff(diff(funcao_green, 1, 1), 1, 2);
    diff_funcao_green(end + 1, :) = diff_funcao_green(end, :);
    diff_funcao_green(:, end + 1) = diff_funcao_green(:, end);

    % Calculando pressao de fato
    superficie_contorno_pressao = abs(Fi_frequencia_escolhida).*abs(diff_funcao_green);
    pressure = 0;
    for y = 1:tamanho_cada_linha
        for x = 1:tamanho_cada_linha
            if y == y1 && x <= x2 || ...
            y <= y2 && x == x2 || ...
            y == y2 && x <= x2 || ...
            y <= y2 && x == x1

                pressure = pressure + superficie_contorno_pressao(y, x);

            end
        end
    end
    pressoes(ponto) = -pressure;
end
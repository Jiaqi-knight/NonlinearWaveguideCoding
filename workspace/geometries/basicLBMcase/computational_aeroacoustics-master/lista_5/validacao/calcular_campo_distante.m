%% calcular_campo_distante: para cada time-step eh um arquivo na pasta 'dados'
function pressoes = calcular_campo_distante(ponto_1_superficie, ponto_2_superficie, ...
U1, U2, frequencia, Angle, raio)
cs = 1/sqrt(3);
cs2 = cs^2;

% fazer o load aqui do tempo total
load('dados.mat');
tempo_total = total_time;

% Efetuando os calculos de F_i
F_i_1 = {};
F_i_2 = {};
for ta = 1:tempo_total
    if  mod(ta, 1000) == 0
        disp('Progresso para os Fi: ');
        disp((ta*100)/tempo_total);
    end

    % Montando os planos
    % invertendo rho linha
    rho = rho_save{ta};
    ux = ux_save{ta};
    uy = uy_save{ta};

    % definindo velocidades
    u1 = ux;
    u2 = uy;
    rho_l = 1;

    Fi_1 = ux;
    Fi_1(:) = 0;
    Fi_2 = ux;
    Fi_2(:) = 0;
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
                Fi_1(posicao_y, posicao_x) = F1;
                Fi_2(posicao_y, posicao_x) = F2;

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
                Fi_1(posicao_y, posicao_x) = F1;
                Fi_2(posicao_y, posicao_x) = F2;

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
                Fi_1(posicao_y, posicao_x) = F1;
                Fi_2(posicao_y, posicao_x) = F2;

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
                Fi_1(posicao_y, posicao_x) = F1;
                Fi_2(posicao_y, posicao_x) = F2;
            end
        end
    end

    %subplot(2,2,1);
    %imagesc(flip((rho - 1)*cs2), [-0.00001 0.00001]);
    %subplot(2,2,2);
    %imagesc(ux, [-0.00001 0.00001]);
    %subplot(2,2,3);
    %imagesc(uy, [-0.00001 0.00001]);
    %subplot(2,2,4);
    %imagesc(Fi, [-0.00001 0.00001]);
    %pause(.0001);
    F_i_1{ta} = Fi_1;
    F_i_2{ta} = Fi_2;
    %size(F_i);
end

% Como ja tenho cada plano ao longo do tempo, agora eh fazer a FFT
%% Calculando Fi_fft
Fi_fft_1 = {};
Fi_fft_2 = {};
for ponto_x_plano = 1:tamanho_cada_linha
    for ponto_y_plano = 1:tamanho_cada_linha
        if ponto_y_plano == y1 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x2 || ...
            ponto_y_plano == y2 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x1

            % Mergulhando nas matrizes
            vetor_auxiliar_1(1:tempo_total) = 0;
            vetor_auxiliar_2(1:tempo_total) = 0;
            for incremento_tempo = 1:tempo_total
                matriz_Fi_1 = F_i_1{incremento_tempo};
                vetor_auxiliar_1(incremento_tempo) = ... 
                matriz_Fi_1(ponto_y_plano, ponto_x_plano);

                matriz_Fi_2 = F_i_2{incremento_tempo};
                vetor_auxiliar_2(incremento_tempo) = ... 
                matriz_Fi_2(ponto_y_plano, ponto_x_plano);
            end

            vetor_auxiliar_1 = vetor_auxiliar_1 - mean(vetor_auxiliar_1);
            vetor_auxiliar_2 = vetor_auxiliar_2 - mean(vetor_auxiliar_2);
            fft_linha_mergulhada_1 = fft(vetor_auxiliar_1);
            fft_linha_mergulhada_2 = fft(vetor_auxiliar_2);
            %figure(1);
            %frequencias = (0:tempo_total-1)/tempo_total;
            %plot(frequencias, abs(fft_linha_mergulhada));
            %figure(2);
            %plot(abs(fft_linha_mergulhada));
            %pause(0.1);

            % Guardando os vetores fft's
            Fi_fft_1{ponto_y_plano, ponto_x_plano} = fft_linha_mergulhada_1/total_time;
            Fi_fft_2{ponto_y_plano, ponto_x_plano} = fft_linha_mergulhada_2/total_time;

        else

        end
    end
end

% Desacoplando os vetores de Fi_fft para montar as superficies certinhas
fft_Fi_1 = {};
fft_Fi_2 = {};
for incremento_tempo = 1:tempo_total
    matriz_auxiliar = F_i_1{incremento_tempo};
    matriz_auxiliar(:) = 0;
    for ponto_y_plano = 1:tamanho_cada_linha
        for ponto_x_plano = 1:tamanho_cada_linha
            if ponto_y_plano == y1 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x2 || ...
            ponto_y_plano == y2 && ponto_x_plano <= x2 || ...
            ponto_y_plano <= y2 && ponto_x_plano == x1

            vetor_fft_auxiliar_1 = Fi_fft_1{ponto_y_plano, ponto_x_plano};
            vetor_fft_auxiliar_2 = Fi_fft_2{ponto_y_plano, ponto_x_plano};
            matriz_auxiliar_1(ponto_y_plano, ponto_x_plano) = ...
            vetor_fft_auxiliar_1(incremento_tempo);
            matriz_auxiliar_2(ponto_y_plano, ponto_x_plano) = ...
            vetor_fft_auxiliar_2(incremento_tempo);

            end
        end
    end
    fft_Fi_1{incremento_tempo} = matriz_auxiliar_1;
    fft_Fi_2{incremento_tempo} = matriz_auxiliar_2;
end

% Selecionando o valor da frequencia
slot_frequencia = 41;
Fi_frequencia_escolhida_1 = fft_Fi_1{slot_frequencia};
Fi_frequencia_escolhida_2 = fft_Fi_2{slot_frequencia};
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
    A = F_i_1{1};
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
    % em y
    diff_funcao_green_y = diff(funcao_green, 1, 1);
    diff_funcao_green_y(end + 1, :) = diff_funcao_green_y(end, :);

    % em x
    diff_funcao_green_x = diff(funcao_green, 1, 2);
    diff_funcao_green_x(:, end + 1) = diff_funcao_green_x(:, end);

    % Calculando pressao de fato
    pressure = 0;
    for y = 1:tamanho_cada_linha
        for x = 1:tamanho_cada_linha
            if y == y1 && x <= x2 || ...
            y <= y2 && x == x2 || ...
            y == y2 && x <= x2 || ...
            y <= y2 && x == x1
                eixo_1 = Fi_frequencia_escolhida_1(y, x)*diff_funcao_green_x(y, x);
                eixo_2 = Fi_frequencia_escolhida_2(y, x)*diff_funcao_green_y(y, x);
                %eixo_1 = trapz(F1*G1 + F2*G2);
                pressure = pressure + abs(eixo_1) + abs(eixo_2);

            end
        end
    end
    pressoes(ponto) = -pressure;
end
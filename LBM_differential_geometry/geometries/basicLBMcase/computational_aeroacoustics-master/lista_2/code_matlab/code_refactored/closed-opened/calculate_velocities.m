function lattice = calculate_velocities(lattice, rho)

    % Determining lattice_time_stephe velocities according to Eq.() (see slides)
    f = lattice{1};
    C_x=[1 0 -1  0 1 -1 -1  1 0];                       % velocity vectors in x
    C_y=[0 1  0 -1 1  1 -1 -1 0];                       % velocity vectors in y
    lattice{2} = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    lattice{3} = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;

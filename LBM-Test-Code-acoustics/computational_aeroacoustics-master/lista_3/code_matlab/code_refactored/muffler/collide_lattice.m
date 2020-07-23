function lattice = collide_lattice(lattice, frequency_relaxation, ... 
    sigma_mat9_left, Ft_left, ...
    sigma_mat9_right, Ft_right);

    density = sum(lattice{1},3);
    w0=16/36. ; w1=4/36. ; w2=1/36.;                    % lattice weights
    rt0= w0*density;
    rt1= w1*density;
    rt2= w2*density;
    velocity_x_pow_2 = lattice{2}.^2; 
    velocity_y_pow_2 = lattice{3}.^2; 
    velocity_pow_2 = velocity_x_pow_2 + velocity_y_pow_2; 
    f1=3.;
    f2=4.5;
    f3=1.5;                                             % coef. of the f equil.
    feq(:,:,1)= rt1 .*(1 +f1*lattice{2} +f2.*velocity_x_pow_2 -f3*velocity_pow_2);
    feq(:,:,2)= rt1 .*(1 +f1*lattice{3} +f2*velocity_y_pow_2 -f3*velocity_pow_2);
    feq(:,:,3)= rt1 .*(1 -f1*lattice{2} +f2*velocity_x_pow_2 -f3*velocity_pow_2);
    feq(:,:,4)= rt1 .*(1 -f1*lattice{3} +f2*velocity_y_pow_2 -f3*velocity_pow_2);
    feq(:,:,5)= rt2 .*(1 +f1*(+lattice{2}+lattice{3}) +f2*(+lattice{2}+lattice{3}).^2 -f3.*velocity_pow_2);
    feq(:,:,6)= rt2 .*(1 +f1*(-lattice{2}+lattice{3}) +f2*(-lattice{2}+lattice{3}).^2 -f3.*velocity_pow_2);
    feq(:,:,7)= rt2 .*(1 +f1*(-lattice{2}-lattice{3}) +f2*(-lattice{2}-lattice{3}).^2 -f3.*velocity_pow_2);
    feq(:,:,8)= rt2 .*(1 +f1*(+lattice{2}-lattice{3}) +f2*(+lattice{2}-lattice{3}).^2 -f3.*velocity_pow_2);
    feq(:,:,9)= rt0 .*(1 - f3*velocity_pow_2);

    f=lattice{1};
    lattice{1} = frequency_relaxation*feq +(1-frequency_relaxation)*f - sigma_mat9_right.*(feq- Ft_right) - sigma_mat9_left.*(feq - Ft_left);
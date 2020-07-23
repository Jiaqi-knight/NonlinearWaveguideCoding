function hydrogenic_orbitals_demo

N = 200;
quantum_number_N = 5;
quantum_number_Ns = {'1','2','3','4','5','6','7','8','9','10'};
quantum_number_L = 'D';
quantum_numbers_Ls = LgivenN(quantum_number_N);
quantum_number_M = -2;
quantum_numbers_Ms = MgivenL(quantum_number_L);
Z = 1;
A = 1;
[w,E_eV,r_mean,r2_mean,a0,a] = Hradial(0,...
    quantum_number_N,...
    orbital2L(quantum_number_L),...
    Z,A);

figure('name','hyrogenic_orbitals_demo','color',[1 1 1] )
surface_handle = make_hydrogenic_orbitals( orbital2L(quantum_number_L), quantum_number_M, N );
title(['N = ',num2str(quantum_number_N),...
    ', L = ',quantum_number_L,...
    ', M = ',num2str(quantum_number_M),...
    ', E = ',num2str(E_eV),'eV'])

%End of code
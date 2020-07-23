%spheria_defaults
% Defines default parameters for the klein spheria. If S.SPHERIA does
% not have an spheria field then it will obtain them from this function.

function d = spheria_defaults

d.N = 200;
d.quantum_number_N = 5;
d.quantum_number_Ns = {'1','2','3','4','5','6','7','8','9','10'};
d.quantum_number_L = 'D';
d.quantum_number_Ls = LgivenN(d.quantum_number_N);
d.quantum_number_M = -2;
d.quantum_number_Ms = MgivenL(d.quantum_number_L);
d.Z = 1;
d.A = 1;
[w,d.E_eV,r_mean,r2_mean,a0,a] = Hradial(0,...
    d.quantum_number_N,...
    orbital2L(d.quantum_number_L),...
    d.Z,d.A);

%Axes properties
d.axes_properties_fields = [];
d.axes_properties = [];

%End of code
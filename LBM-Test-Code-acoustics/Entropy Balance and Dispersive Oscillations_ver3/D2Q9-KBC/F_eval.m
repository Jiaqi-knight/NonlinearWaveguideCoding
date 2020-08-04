function F=F_eval(scheme,populations,popequilibriums,alfa)
%-----------------------------------------------------------
%    Function that evaluate the following function:
%    F(alfa)=H(f_i)-H(f_i+alfa*(f_i-f_i^eq)) Dorschner-2012-17
%-----------------------------------------------------------
H1=Hf(scheme,populations);
H2=Hf(scheme,populations+alfa*(popequilibriums-populations));
F=H1-H2;

end
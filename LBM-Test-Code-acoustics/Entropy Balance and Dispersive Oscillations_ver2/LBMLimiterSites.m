% A function to detect an appropriate lattice site for limiting and to
% measure it's position relative to the leading edge of the shock.
function [Sites,shockRelative] = LBMLimiterSites(populations,...
popequilibriums,nx,Limiter)
shockRelative = 0;
if (Limiter == 0)
Sites = [];
elseif (Limiter == 1)
Sf = LBMEvaluateS(populations,popequilibriums,zeros(1,nx));
Sfeq = LBMEvaluateS(popequilibriums,popequilibriums,zeros(1,nx));
DeltaS = Sfeq - Sf;
shockPos = max(find(DeltaS > 10^-15));
[M,I] = max(DeltaS);
Sites = I;
shockRelative = shockPos - Sites;
end

% A function to find an interior approximation to the root of the entropy
% parabola
function intPab = LBMInteriorParabola(scheme,populations,popequilibriums, ...
STarget,alpha)
SAlphaZero = LBMEvaluateS(scheme,populations,popequilibriums,alpha);
SDashAlphaZero = LBMEvaluateDiffS(scheme,populations,popequilibriums,alpha);
SDash2AlphaZero = LBMEvaluateDiff2S(scheme,populations,popequilibriums,alpha);                    
%eq-3.103
% Pab_min=max(populations./(populations-popequilibriums))
% Pab_min=

intPab = [];
for j = 1:length(alpha)
intPab(j) = max(roots([0.5.*SDash2AlphaZero(j) SDashAlphaZero(j) ...
SAlphaZero(j) - STarget(j)])) + alpha(j);
end
if(~isreal(intPab))
%    keyboard
   intPab=real(intPab);
end
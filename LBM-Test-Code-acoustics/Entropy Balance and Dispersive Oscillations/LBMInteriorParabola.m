% A function to find an interior approximation to the root of the entropy
% parabola
function intPab = LBMInteriorParabola(populations,popequilibriums, ...
STarget,alpha)
SAlphaZero = LBMEvaluateS(populations,popequilibriums,alpha);
SDashAlphaZero = LBMEvaluateDiffS(populations,popequilibriums,alpha);
SDash2AlphaZero = LBMEvaluateDiff2S(populations,popequilibriums,alpha);                    
intPab = [];
for j = 1:length(alpha)
intPab(j) = max(roots([0.5.*SDash2AlphaZero(j) SDashAlphaZero(j) ...
SAlphaZero(j) - STarget(j)])) + alpha(j);
end
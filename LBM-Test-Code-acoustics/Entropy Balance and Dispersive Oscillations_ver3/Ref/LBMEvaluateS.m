% A function to evaluate the entropy.
function S = LBMEvaluateS(scheme,populations,popequilibriums,alpha)
alphapop = populations + (ones(size(populations,1),1)*alpha).*(popequilibriums ...
- populations);
w=scheme(:,2);
% S = -alphapop(1,:).*log(alphapop(1,:)) ...
% - alphapop(2,:).*log(alphapop(2,:)./4) ...
% - alphapop(3,:).*log(alphapop(3,:)); %eq-2.5

S = - sum(alphapop.*log(alphapop./w),1);
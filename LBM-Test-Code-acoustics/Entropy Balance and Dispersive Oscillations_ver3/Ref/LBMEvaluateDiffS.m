% A function to find the first derivative of the entropy
function DiffS = LBMEvaluateDiffS(scheme,populations,popequilibriums,alpha)
alphapop = populations + (ones(size(scheme,1),1)*alpha).*(popequilibriums ...
- populations);
partDiff = popequilibriums - populations;
% DiffS = -partDiff(1,:).*(log(alphapop(1,:)) + 1) ...
% -partDiff(2,:).*(log(alphapop(2,:)./4) + 1) ...
% -partDiff(3,:).*(log(alphapop(3,:)) + 1);

w=scheme(:,2);
DiffS = - sum(partDiff.*(log(alphapop./w)+1),1);
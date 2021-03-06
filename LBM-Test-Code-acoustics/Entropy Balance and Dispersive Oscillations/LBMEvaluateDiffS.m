% A function to find the first derivative of the entropy
function DiffS = LBMEvaluateDiffS(populations,popequilibriums,alpha)
alphapop = populations + (ones(3,1)*alpha).*(popequilibriums ...
- populations);
partDiff = popequilibriums - populations;
DiffS = -partDiff(1,:).*(log(alphapop(1,:)) + 1) ...
-partDiff(2,:).*(log(alphapop(2,:)./4) + 1) ...
-partDiff(3,:).*(log(alphapop(3,:)) + 1);
% A function to find the second derivative of the entropy.
function Diff2S = LBMEvaluateDiff2S(populations,popequilibriums,alpha)
alphapop = populations + (ones(3,1)*alpha).*(popequilibriums - populations);
partDiff = popequilibriums - populations;
Diff2S = - partDiff(1,:).^2./(alphapop(1,:)) ...
- partDiff(2,:).^2./(alphapop(2,:)) ...
- partDiff(3,:).^2./(alphapop(3,:));

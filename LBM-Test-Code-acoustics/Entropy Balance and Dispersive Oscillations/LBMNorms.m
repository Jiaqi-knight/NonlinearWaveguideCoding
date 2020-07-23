% A function to implement the norms necessary to measure convergence of the
% root finding
function norms = LBMNorms(populations,popequilibriums,nx,nChoice);
norms = zeros(1,nx);
if (nChoice == 0)
for j = 1:nx
norms(j) = norm(populations(:,j) - popequilibriums(:,j),1);
end
elseif(nChoice == 1)
for j = 1:nx
poptemp = (populations(1,j) ...
- popequilibriums(1,j)).^2./popequilibriums(1,j);
poptemp = poptemp + (populations(2,j) ...
- popequilibriums(2,j)).^2./popequilibriums(2,j);
poptemp = poptemp + (populations(3,j) ...
- popequilibriums(3,j)).^2./popequilibriums(3,j);
norms(j) = sqrt(poptemp);
end
end
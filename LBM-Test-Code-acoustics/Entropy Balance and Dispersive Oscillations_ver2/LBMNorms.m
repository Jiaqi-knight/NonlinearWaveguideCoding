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
        poptemp=0;
        for q=1:size(populations,1)
            poptemp=poptemp+(populations(q,j) ...
                - popequilibriums(q,j)).^2./popequilibriums(q,j);
            norms(j) = sqrt(poptemp);
        end
    end
end
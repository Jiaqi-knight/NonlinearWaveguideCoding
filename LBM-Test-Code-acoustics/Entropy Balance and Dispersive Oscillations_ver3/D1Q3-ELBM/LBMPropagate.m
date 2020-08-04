%A function to propagate the populations
function populations = LBMPropagate(oldpopulations,nx,scheme)
populations = oldpopulations;
c=scheme(:,1);
for k=1:length(c)
    if c(k)>0
%         left=populations(k,nx-c(k)+1:nx); %periodic condition
        populations(k,c(k)+1:nx)   = populations(k,1:nx-c(k));
%         populations(k,1:c(k))  =  left;
    elseif c(k)<0
%         right=populations(k,1:-c(k));
        populations(k,1:nx+c(k)) = populations(k,1-c(k):nx);
%         populations(k,nx+c(k)+1:nx) = right;
    end
end
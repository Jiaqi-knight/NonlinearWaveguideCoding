%A function to propagate the populations
function populations = LBMPropagate(oldpopulations,nx,ny,scheme)
populations = oldpopulations;
c_x=scheme(:,1);
c_y=scheme(:,2);

for k=1:length(c_x)
    if c_x(k)>0
%         left=populations(k,nx-c(k)+1:nx); %periodic condition
        populations(c_x(k)+1:nx,:,k)   = populations(1:nx-c_x(k),:,k);
%         populations(k,1:c(k))  =  left;
    elseif c_x(k)<0
%         right=populations(k,1:-c(k));
        populations(1:nx+c_x(k),:,k) = populations(1-c_x(k):nx,:,k);
%         populations(k,nx+c(k)+1:nx) = right;
    end
end
for k=1:length(c_y)
    if c_y(k)>0
%         left=populations(k,nx-c(k)+1:nx); %periodic condition
        populations(:,c_y(k)+1:ny,k)   = populations(:,1:ny-c_y(k),k);
%         populations(k,1:c(k))  =  left;
    elseif c_y(k)<0
%         right=populations(k,1:-c(k));
        populations(:,1:ny+c_y(k),k) = populations(:,1-c_y(k):ny,k);
%         populations(k,nx+c(k)+1:nx) = right;
    end
end
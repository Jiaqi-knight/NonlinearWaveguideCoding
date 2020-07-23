function [rho, u] = Macroscopic(f,c)
% Compute macroscopic quantities
fu = zeros(size(f,1),1);

% Density
rho = sum(f,2);

% Mass flux
for k = 1:length(c)
    fu = fu + (c(k) * f(:,k));
end


% Get velocity
u = fu ./ rho;


end
function [P, J] = Macroscopic(f,xi,LatTyp)
% Compute macroscopic quantities
fu = zeros(1,size(f,2),size(f,3),size(f,4));
fv = zeros(1,size(f,2),size(f,3),size(f,4));
fw = zeros(1,size(f,2),size(f,3),size(f,4));

% Density
P(:,:,:) = sum(f,1);

% Mass flux
for v = 1:7
    fu = fu + (xi(1,v) * f(v,:,:,:));
    fv = fv + (xi(2,v) * f(v,:,:,:));
    fw = fw + (xi(3,v) * f(v,:,:,:));
end

% % Add forces
% fu = fu + (0.5 * reshape(P,[1 size(P)]) .* force_xyz(1,:,:,:));
% fv = fv + (0.5 * reshape(P,[1 size(P)]) .* force_xyz(2,:,:,:));
% fw = fw + (0.5 * reshape(P,[1 size(P)]) .* force_xyz(3,:,:,:));

% Get velocity
J(1,:,:,:) = fu ;
J(2,:,:,:) = fv ;
J(3,:,:,:) = fw ;

% Apply resets
J( repmat( strcmp(LatTyp,{'b'}), 3,1,1,1 ) ) = 0;

% % Total kinetic energy
% ke = 0.5*sum(sum(sum(sum( sqrt(J(1,:,:,:).^2 + J(2,:,:,:).^2 + J(3,:,:,:).^2) ))));


end
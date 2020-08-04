function [rho,u,f] = initialfield(rho,u)
% This function is to Initial the flow field
% Inputs is the macroscopic field varibles rho(Nx,Ny), u(Nx,Ny,D) and
% microscopic varible f(Nx,Ny,Q)
% Outputs is also the macroscopic field varibles rho(Nx,Ny), u(Nx,Ny,D) and
% microscopic field varibles f(Nx,Ny,Q)

global Q rho_0 U Nx Ny
% Macroscopic
rho(:,:) = rho_0;     % initial the density field 
u(:,:,:) = 0;         % initial the velocity field (Inner)
u(:,Ny+1,1) = U;      % initial the velocity field (Boundary Lid as a constant velocity)

% Mesoscopic
% initial the distribution function  f_initial = f_eq
f = zeros(Nx+1,Ny+1,Q);    % preallocated f(Nx,Ny,Q) 
temp_feq = zeros(Q,1);     % preallocated temporary feq
for i  = 1:Nx+1
    for j = 1:Ny+1
        temp_feq = f_eq(rho(i,j),u(i,j,:));
        f(i,j,:) = temp_feq; 
    end
end


end

% tiem O(Nx*Ny)

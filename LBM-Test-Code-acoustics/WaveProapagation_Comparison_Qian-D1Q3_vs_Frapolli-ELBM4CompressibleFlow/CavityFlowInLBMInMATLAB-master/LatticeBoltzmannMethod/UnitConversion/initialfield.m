function [rho,ux,uy,f] = initialfield()
% This function is to Initial the flow field
% Inputs is the macroscopic field varibles rho(Nx,Ny), u(Nx,Ny,D) and
% microscopic varible f(Nx,Ny,Q)
% Outputs is also the macroscopic field varibles rho(Nx,Ny), u(Nx,Ny,D) and
% microscopic field varibles f(Nx,Ny,Q)

global  Q rho_0 U omega Nx Ny
% Macroscopic
rho = zeros(Nx+1,Ny+1);
rho(:,:) = rho_0;           % initial the density field 
ux = zeros(Nx+1,Ny+1);      % initial the velocity field (Inner)
uy = ux;
ux(:,Ny+1) = U;             % initial the velocity field (Boundary Lid as a constant velocity)

% Mesoscopic
% initial the distribution function  f_initial = f_eq (Using vector operator)
f = reshape(ones( (Nx+1)*(Ny+1) , 1) * omega, Nx+1, Ny+1, Q);
[f(:,Ny+1,1),f(:,Ny+1,2),f(:,Ny+1,3),f(:,Ny+1,4),f(:,Ny+1,5),f(:,Ny+1,6),...
    f(:,Ny+1,7),f(:,Ny+1,8),f(:,Ny+1,9)]= arrayfun(@f_eq,rho(:,Ny+1),ux(:,Ny+1),uy(:,Ny+1));

% initial the distribution function  f_initial = f_eq (Using for loop)
% f = zeros(Nx+1,Ny+1,Q);    % preallocated f(Nx,Ny,Q) 
% temp_feq = zeros(Q,1);     % preallocated temporary feq
% for i  = 1:Nx+1
%     for j = 1:Ny+1
%         temp_feq = f_eq(rho(i,j),u(i,j,:));
%         f(i,j,:) = temp_feq; 
%     end
% end


end

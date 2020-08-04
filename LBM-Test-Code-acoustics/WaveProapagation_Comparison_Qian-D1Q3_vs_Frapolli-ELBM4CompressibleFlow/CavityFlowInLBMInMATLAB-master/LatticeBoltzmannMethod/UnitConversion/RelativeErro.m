function varepsilon = RelativeErro(ux, uy, ux_temp, uy_temp)
% This function is designed to compute the relative erro in velocity field
% of two near time step.
%  inputs field varibles is u(Nx+1,Ny+1,2), u_temp(Nx+1,Ny+1,2)
% to addtion the whole field's erro
num = (ux - ux_temp).^2 + (uy - uy_temp).^2; % Square
num = sum(sum(num));                         % Sum
den = ux_temp.^2 + uy_temp.^2;               % Square
den = sum(sum(den));                         % Sum
varepsilon = sqrt(num)/sqrt(den) ;    
end


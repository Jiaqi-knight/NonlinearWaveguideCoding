function varepsilon = RelativeErro(u, u_temp)
% This function is designed to compute the relative erro in velocity field
% of two near time step.
%  inputs field varibles is u(Nx+1,Ny+1,2), u_temp(Nx+1,Ny+1,2)
global Nx Ny Dimension
varepsilon = 0;
num = 0;
den = 0;
% to addtion the whole field's erro
for i = 1:Nx+1
    for j = 1:Ny+1
        for m = 1:Dimension
            num = num + (u(i,j,m) - u_temp(i,j,m))^2;
            den = den + u_temp(i,j,m)^2;
        end
    end
end
varepsilon = sqrt(num)/sqrt(den) ;    
end


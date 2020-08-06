% Macroscopic quantity calculation
function LBM_macro( r ) 
global Grids dims nVels c
% Declarations
N_lim =length( Grids{r}.XPos );
M_lim =length( Grids{r}.YPos );
K_lim =length( Grids{r}.ZPos );
rho_temp = 0;
fux_temp = 0;
fuy_temp = 0;
fuz_temp = 0;

% Loop over lattice
for ( i = 1 : N_lim )
    for ( j = 1 : M_lim )
        for ( k = 1 : K_lim )
            
            % Reset temporary variables
            rho_temp = 0; fux_temp = 0; fuy_temp = 0; fuz_temp = 0;
            
            % Get index
            idxrho = idxmap(i,j,k,M_lim,K_lim);
            
            for ( v = 1 : nVels)
                
                % Get index
                idxf = idxmap(i,j,k,v,M_lim,K_lim,nVels);
                
                % Sum up to find momentum
                fux_temp =fux_temp+ c(1,v) * Grids{r}.f(idxf);
                fuy_temp =fuy_temp+ c(2,v) * Grids{r}.f(idxf);
                fuz_temp =fuz_temp+ c(3,v) * Grids{r}.f(idxf);
                
                % Sum up to find density
                rho_temp =rho_temp+ Grids{r}.f(idxf);
                
            end
            
            % Assign density
            Grids{r}.rho(idxrho) = rho_temp;
            
            % Get index
            idx_ux = idxmap(i,j,k,1,M_lim,K_lim,dims);
            idx_uy = idxmap(i,j,k,2,M_lim,K_lim,dims);
            idx_uz = idxmap(i,j,k,3,M_lim,K_lim,dims);
            
            % Assign velocity
            Grids{r}.u(idx_ux) = fux_temp / rho_temp;
            Grids{r}.u(idx_uy) = fuy_temp / rho_temp;
            if (dims == 3)
                Grids{r}.u(idx_uz) = fuz_temp / rho_temp;
            end
            
        end
    end
end


end
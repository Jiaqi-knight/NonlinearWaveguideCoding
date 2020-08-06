% Streaming operator
% Applies periodic BCs on level 0
function LBM_stream( r )
global Grids nVels c
% Declarations
N_lim = length(Grids{r}.XPos);
M_lim = length(Grids{r}.YPos);
K_lim = length(Grids{r}.ZPos);
% dest_x, dest_y, dest_z;

% Create temporary lattice of zeros to prevent overwriting useful populations
% vector<double> f_new( Grids[r].f.size(), 0.0 );

% Stream one lattice site at a time
for ( i = 1 : N_lim )
    for ( j = 1 : M_lim )
        for ( k = 1 : K_lim )
            
            for ( v = 1 : nVels )
                
                % Get index
                idx_ijk = idxmap(i,j,k,M_lim,K_lim);
                
                % If fine site then do not stream in any direction
                if (Grids{r}.LatTyp(idx_ijk) == 2)
                    break;
                end
                
                % Only apply periodic BCs on coarsest level
                if (r == 1)
                    % Compute destination coordinates
                    dest_x = mod((i+c(1,v) + N_lim-1) , N_lim)+1;
                    dest_y = mod((j+c(2,v) + M_lim-1) , M_lim)+1;
                    dest_z = mod((k+c(3,v) + K_lim-1) , K_lim)+1;
                else
                    % No periodic BCs
                    dest_x = i+c(1,v);
                    dest_y = j+c(2,v);
                    dest_z = k+c(3,v);
                end
                
                % Get destination index
                idx_dest = idxmap(dest_x,dest_y,dest_z,M_lim,K_lim);
                
                % If destination off-grid, do not stream
                if (	(dest_x >= N_lim+1 || dest_x < 1) ||...
                        (dest_y >= M_lim+1 || dest_y < 1) ||...
                        (dest_z >= K_lim+1 || dest_z < 1)...
                        )
                    % Do nothing
                    
                else
                    
                    % Check destination site type and decide whether to stream or not
                    if ( (Grids{r}.LatTyp(idx_dest) == 2) ||...
                            ( (Grids{r}.LatTyp(idx_dest) == 4)...
                            && (Grids{r}.LatTyp(idx_ijk) == 4) ) )
                        % Fine -- ignore
                        % TL lower level to TL lower level -- done on lower grid stream so ignore too.
                        
                    else
                        
                        % Get index
                        idx_destv = idxmap(dest_x,dest_y,dest_z,v,M_lim,K_lim,nVels);
                        idx_ijkv = idxmap(i,j,k,v,M_lim,K_lim,nVels);
                        
                        % Stream population
                        f_new(idx_destv,1) = Grids{r}.f(idx_ijkv);
                        
                    end
                end
                
            end
            
        end
    end
end

% Replace old grid with new grid
Grids{r}.f = f_new;

end

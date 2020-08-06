% Explosion operation -- called from fine level
function LBM_explode(  r ) 
global nVels Grids dims
% Declarations
% y_start, x_start, z_start;
M_fine =length( Grids{r}.YPos);
K_fine =length( Grids{r}.ZPos);
M_coarse = length(Grids{r-1}.YPos);
K_coarse = length(Grids{r-1}.ZPos);

% Loop over coarse grid
for (i = 1 : length(Grids{r-1}.XPos))
    for (j = 1 : length(Grids{r-1}.YPos))
        for ( k = 1 : length(Grids{r-1}.ZPos))
            
            idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);
            
            % If TL to lower level then partitioning required
            if (Grids{r-1}.LatTyp(idx_coarse) == 4)
                
                % Lookup indices for lower level. L0 start is different from other
                % levels as it only has a single TL.
                if (r-1 == 1) 
                    x_start = Grids{1}.XInd(1);
                    y_start = Grids{1}.YInd(1);
                    z_start = Grids{1}.ZInd(1);
                     else 
                        x_start = Grids{r-1}.XInd(3);
                        y_start = Grids{r-1}.YInd(3);
                        if (dims == 3)
                        z_start = Grids{r-1}.ZInd(3);
                        else
                            z_start = Grids{r-1}.ZInd(1); % if 2D set to default
                        end
                end
                
                % Update fine grid values according to Rohde et al.
                for ( v = 1 : nVels)
                    
                    % Get coarse site value
                    idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);
                    coarse_f = Grids{r-1}.f(idx_coarsef);
                    
                    % Find indices of fine site
                    idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
                    fi = idx_fine(1);
                    fj = idx_fine(2);
                    fk = idx_fine(3);
                    
                    if (dims == 3)
                        
                        % 3D Case -- cube of 8 cells
                        
                        % Flatten indices for each cell in the cube
                        idx1 = idxmap(fi,	fj,		fk,		v,M_fine,K_fine,nVels);
                        idx2 = idxmap(fi+1,	fj,		fk,		v,M_fine,K_fine,nVels);
                        idx3 = idxmap(fi,	fj+1,	fk,		v,M_fine,K_fine,nVels);
                        idx4 = idxmap(fi+1,	fj+1,	fk,		v,M_fine,K_fine,nVels);
                        idx5 = idxmap(fi,	fj,		fk+1,	v,M_fine,K_fine,nVels);
                        idx6 = idxmap(fi+1,	fj,		fk+1,	v,M_fine,K_fine,nVels);
                        idx7 = idxmap(fi,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
                        idx8 = idxmap(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
                        
                        % Copy coarse to fine
                        Grids{r}.f(idx1) = coarse_f;
                        Grids{r}.f(idx2) = coarse_f;
                        Grids{r}.f(idx3) = coarse_f;
                        Grids{r}.f(idx4) = coarse_f;
                        Grids{r}.f(idx5) = coarse_f;
                        Grids{r}.f(idx6) = coarse_f;
                        Grids{r}.f(idx7) = coarse_f;
                        Grids{r}.f(idx8) = coarse_f;
                        
                    else
                        
                        % 2D Case -- square of 4 cells
                        
                        % Flatten indices
                        idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
                        idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
                        idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
                        idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
                        
                        % Copy coarse to fine
                        Grids{r}.f(idx1) = coarse_f;
                        Grids{r}.f(idx2) = coarse_f;
                        Grids{r}.f(idx3) = coarse_f;
                        Grids{r}.f(idx4) = coarse_f;
                        
                    end
                end
                
            end
            
        end
    end
end


end


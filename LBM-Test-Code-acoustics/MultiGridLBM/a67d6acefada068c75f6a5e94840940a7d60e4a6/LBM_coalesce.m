% Coalesce operation -- called from coarse level
function LBM_coalesce(  r  )
global Grids dims nVels
% Declarations
% 	 y_start, x_start, z_start;
M_coarse = length(Grids{r}.YPos);
K_coarse = length(Grids{r}.ZPos);
M_fine = length(Grids{r+1}.YPos);
K_fine = length(Grids{r+1}.ZPos);

% Loop over coarse grid
for (i = 1:length(Grids{r}.XPos) )
    for (j = 1 : length(Grids{r}.YPos) )
        for (k = 1 : length(Grids{r}.ZPos) )
            
            idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);
            
            % If TL to lower level then fetch values from lower level
            if (Grids{r}.LatTyp(idx_coarse) == 4)
                
                % Lookup indices for lower level. L0 start is different from other
                % levels as it only has a single TL.
                if (r == 1)
                    x_start = Grids{1}.XInd(1);
                    y_start = Grids{1}.YInd(1);
                    z_start = Grids{1}.ZInd(1);
                else
                    x_start = Grids{r}.XInd(3);
                    y_start = Grids{r}.YInd(3);
                    if (dims == 3)
                        z_start = Grids{r}.ZInd(3);
                    else
                        z_start = Grids{r}.ZInd(1); % if 2D set to default
                    end
                end
                
                % Loop over directions
                for ( v = 1 : nVels )
                    
                    % Get coarse site index
                    idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);
                    
                    % Check to see if f value is missing on coarse level
                    if (Grids{r}.f(idx_coarsef) ~= 0)
                        % If not do nothing
                        
                    else
                        
                        % Find indices of corresponding fine site
                        idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
                        fi = idx_fine(1);
                        fj = idx_fine(2);
                        fk = idx_fine(3);
                        
                        if (dims == 3)
                            
                            % 3D Case -- cube of 8 cells
                            
                            % Now flatten indices for each cell in the cube
                            idx1 = idxmap(fi,	fj,		fk,		v,M_fine,K_fine,nVels);
                            idx2 = idxmap(fi+1,	fj,		fk,		v,M_fine,K_fine,nVels);
                            idx3 = idxmap(fi,	fj+1,	fk,		v,M_fine,K_fine,nVels);
                            idx4 = idxmap(fi+1,	fj+1,	fk,		v,M_fine,K_fine,nVels);
                            idx5 = idxmap(fi,	fj,		fk+1,	v,M_fine,K_fine,nVels);
                            idx6 = idxmap(fi+1,	fj,		fk+1,	v,M_fine,K_fine,nVels);
                            idx7 = idxmap(fi,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
                            idx8 = idxmap(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
                            
                            % Average the values
                            Grids{r}.f(idx_coarsef) = (...
                                Grids{r+1}.f(idx1) + Grids{r+1}.f(idx2) + Grids{r+1}.f(idx3) + Grids{r+1}.f(idx4) ...
                                + Grids{r+1}.f(idx5) + Grids{r+1}.f(idx6) + Grids{r+1}.f(idx7) + Grids{r+1}.f(idx8)...
                                ) / power(2, dims);
                            
                        else
                            
                            % 2D Case -- square of 4 cells
                            
                            % Flatten indices
                            idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
                            idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
                            idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
                            idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
                            
                            % Average the values
                            Grids{r}.f(idx_coarsef) = (...
                                Grids{r+1}.f(idx1) + Grids{r+1}.f(idx2) + Grids{r+1}.f(idx3) + Grids{r+1}.f(idx4) ...
                                ) / power(2, dims);
                            
                        end
                    end
                    
                end
                
            end
            
        end
    end
end

end


% ***************************************************************************************************
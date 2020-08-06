%% Generate all the quantities for the extra (refined) grids
%% Assumes a volumetric setup
function LBM_init_multi (  )

%% Generate NODE NUMBERS
%% Check grid core to make sure it can support a finer grid, if not throw exception
global Grids dims Nref N M K nVels


if (dims == 3)
    
    if ( (length(Grids{1}.XInd) < 3 || length(Grids{1}.YInd) < 3 || length(Grids{1}.ZInd) < 3) ||...
            ((length(Grids{1}.XInd) == 3 || length(Grids{1}.YInd) == 3 || length(Grids{1}.ZInd) == 3)...
            && Nref > 1) )
        disp('error');
    else
        %% First level unique as doesn't need TL
        Grids{2}.XInd = [2: length(Grids{1}.XInd)*2 ];
        Grids{2}.YInd = [2: length(Grids{1}.YInd)*2 ];
        Grids{2}.ZInd = [2: length(Grids{1}.ZInd)*2 ];
    end
    
elseif (dims == 2)
    
    if ( (length(Grids{1}.XInd) < 3 || length(Grids{1}.YInd) < 3) ||...
            ((length(Grids{1}.XInd) == 3 || length(Grids{1}.YInd) == 3)...
            && Nref > 1) )
        disp('error');
    else
        %% First level unique as doesn't need TL
        Grids{2}.XInd = [1: length(Grids{1}.XInd)*2 ];
        Grids{2}.YInd = [1: length(Grids{1}.YInd)*2 ];
        Grids{2}.ZInd = 1; %% Default for 2D
        
    end
    
end

%% Lower levels
if (Nref > 1)
    for ( r = 2 : Nref)
        %% Increase refinement by factor of 2
        Grids{r}.XInd = [1: (length(Grids{r}.XInd)-4)*2 ];
        Grids{r}.YInd = [1: (length(Grids{r}.YInd)-4)*2 ];
        if (dims == 3)
            Grids{r}.ZInd = [1: (length(Grids{r}.ZInd)-4)*2 ];
        elseif (dims == 2)
            Grids{r}.ZInd = 1; %% Default for 2D
        end
    end
end





%% Correct L0 TYPING MATRIX to add suitable labels
if (dims == 3)
    for (i = Grids{1}.XInd(1): Grids{1}.XInd(end-1)+1)
        for ( j = Grids{1}.YInd(1): Grids{1}.YInd(end-1)+1)
            for ( k = Grids{1}.ZInd(1): Grids{1}.ZInd(end-1)+1)
                
                idx = idxmap(i,j,k,M,K);
                Grids{1}.LatTyp(idx) = 4;
                
            end
        end
    end
    for (i = Grids{1}.XInd(2): Grids{1}.XInd(end-2)+1)
        for ( j = Grids{1}.YInd(2): Grids{1}.YInd(end-2)+1)
            for ( k = Grids{1}.ZInd(2): Grids{1}.ZInd(end-2)+1)
                
                idx = idxmap(i,j,k,M,K);
                Grids{1}.LatTyp(idx) = 2;
            end
        end
    end
    
else
    for ( i = Grids{1}.XInd(1): Grids{1}.XInd(end-1)+1)
        for ( j = Grids{1}.YInd(1):  Grids{1}.YInd(end-1)+1)
            k = 1;
            
            idx = idxmap(i,j,k,M,K);
            Grids{1}.LatTyp(idx) = 4;
            
        end
    end
    for ( i = Grids{1}.XInd(2): Grids{1}.XInd(end-2)+1)
        for ( j = Grids{1}.YInd(2):  Grids{1}.YInd(end-2)+1)
            k = 1;
            
            idx = idxmap(i,j,k,M,K);
            Grids{1}.LatTyp(idx) = 2;
            
        end
    end
end

%% Generate lower level TYPING MATRICES
if (dims == 3)
    for ( r = 1: Nref)
        
        %% Get grid sizes
        M_lim = length(Grids{r+1}.YInd);
        K_lim = length(Grids{r+1}.ZInd);
        
        %% Resize
        Grids{r+1}.LatTyp=zeros(length(Grids{r+1}.YInd)*length(Grids{r+1}.XInd)*length(Grids{r+1}.ZInd),1);
        
        %% Start with TL from level above
        for (i = 1: length(Grids{r+1}.XInd))
            for (j = 1: length(Grids{r+1}.YInd))
                for (k = 1: length(Grids{r+1}.ZInd))
                    
                    idx = idxmap(i,j,k,M_lim,K_lim);
                    Grids{r+1}.LatTyp(idx) = 3;
                    
                end
            end
        end
        
        %% Check if lower grids exist
        if (Nref > r)
            
            %% Add TL for next level down
            for (i = 3 :  length(Grids{r}.XInd)-2)
                for (j = 3 : length(Grids{r}.YInd)-2)
                    for (k = 3 : Grids{r}.ZInd.size()-2)
                        
                        idx = idxmap(i,j,k,M_lim,K_lim);
                        Grids{r}.LatTyp(idx) = 4;
                        
                    end
                end
            end
            
            
            %% Label rest as fine
            for (i = 4: length(Grids{r}.XInd)-3)
                for (j = 4:  length(Grids{r}.YInd)-3)
                    for (k = 4: Grids{r}.ZInd.size()-3)
                        
                        idx = idxmap(i,j,k,M_lim,K_lim);
                        Grids{r}.LatTyp(idx) = 2;
                        
                    end
                end
            end
            
            
        else
            %% Reached lowest level so label rest as coarse
            for (i = 3: length(Grids{r}.XInd)-2)
                for (j = 3:  length(Grids{r}.YInd)-2)
                    for (k = 3:  Grids{r}.ZInd.size()-2)
                        
                        idx = idxmap(i,j,k,M_lim,K_lim);
                        Grids{r+1}.LatTyp(idx) = 1;
                        
                    end
                end
            end
            
        end
    end
    
    
else
    for ( r = 1: Nref)
        
        %% Get grid sizes
        M_lim = length(Grids{r+1}.YInd);
        K_lim = length(Grids{r+1}.ZInd);
        
        %% Resize
        Grids{r+1}.LatTyp=zeros(length(Grids{r+1}.YInd)*length(Grids{r+1}.XInd)*length(Grids{r+1}.ZInd),1);
        
        %% Start with TL from level above
        for (i = 1 : length(Grids{r+1}.XInd))
            for (j = 1 : length(Grids{r+1}.YInd))
                k = 1;
                
                idx = idxmap(i,j,k,M_lim,K_lim);
                Grids{r+1}.LatTyp(idx) = 3;
                
            end
        end
        
        %% Check if lower grids exist
        if (Nref > r)
            
            %% Add TL for next level down
            for (i = 3 : length(Grids{r+1}.XInd)-2)
                for (j = 3 : length(Grids{r+1}.YInd)-2)
                    k = 1;
                    
                    idx = idxmap(i,j,k,M_lim,K_lim);
                    Grids{r+1}.LatTyp(idx) = 4;
                    
                end
            end
            
            
            %% Label rest as fine
            for (i= 4 : length(Grids{r}.XInd)-3)
                for (j= 4 : length(Grids{r}.YInd)-3)
                    k = 1;
                    
                    idx = idxmap(i,j,k,M_lim,K_lim);
                    Grids{r}.LatTyp(idx) = 2;
                    
                end
            end
            
        else
            %% Reached lowest level so label rest as coarse
            for (i = 3 : length(Grids{r+1}.XInd)-2)
                for (j = 3 : length(Grids{r+1}.YInd)-2)
                    k = 1;
                    
                    idx = idxmap(i,j,k,M_lim,K_lim);
                    Grids{r+1}.LatTyp(idx) = 1;
                    
                end
            end
            
            
        end
    end
    
    %% Generate POSITION VECTORS of nodes
%     count, first_idx, last_idx;
    
    %% Define spacing
    Grids{2}.dx = Grids{1}.dx/2;
    Grids{2}.dy = Grids{1}.dy/2;
    Grids{2}.dz = Grids{1}.dz/2;
    
    %% L0 unusual due to no extra TL
    count = 0;
    first_idx = Grids{1}.XInd(1);
    last_idx = Grids{1}.XInd( length(Grids{1}.XInd) );
    for ( i = first_idx :  last_idx )
        
        %% Call position mapping function
        val1 = posmapref(Grids{1}.XPos(i),2,'x','+');
        val2 = posmapref(Grids{1}.XPos(i),2,'x','-');
        
        %% Assign to Position vectors
        Grids{2}.XPos(1 + count) = val2 ;
        Grids{2}.XPos(1 + count+1) = val1 ;
        
        %% Increment counter
        count=count+2;
        
    end
    
    count = 0;
    first_idx = Grids{1}.YInd(1);
    last_idx = Grids{1}.YInd( length(Grids{1}.YInd) );
    for ( j = first_idx : last_idx )
        
        %% Call position mapping function
        val1 = posmapref(Grids{1}.YPos(j),2,'y','+');
        val2 = posmapref(Grids{1}.YPos(j),2,'y','-');
        
        %% Assign to Position vectors
        Grids{2}.YPos(1 + count) = val2 ;
        Grids{2}.YPos(1 + count+1) = val1 ;
        
        %% Increment counter
        count =count+2;
        
    end
    
    if (dims == 3)
        count = 0;
        first_idx = Grids{1}.ZInd(1);
        last_idx = Grids{1}.ZInd( length(Grids{1}.ZInd) );
        for ( k = first_idx : last_idx)
            
            %% Call position mapping function
            val1 = posmapref(Grids{1}.ZPos(k),2,'z','+');
            val2 = posmapref(Grids{1}.ZPos(k),2,'z','-');
            
            %% Assign to Position vectors
            Grids{2}.ZPos(1 + count) = val2 ;
            Grids{2}.ZPos(1 + count+1) = val1 ;
            %% Increment counter
            count=count+2;
            
        end
        
    else
        Grids{2}.ZPos=1; %% 2D default
    end
    
    
    %% Now do lower levels
    for ( r = 2 : Nref)
        
        %% Spacing
        Grids{r+1}.dx = Grids{r}.dx/2;
        Grids{r+1}.dy = Grids{r}.dy/2;
        Grids{r+1}.dz = Grids{r}.dz/2;
        
        count = 0;
        first_idx = Grids{r}.Grids(3);
        last_idx = Grids{r}.XInd( length(Grids{r}.XInd)-2 );
        for ( i = first_idx : last_idx)
            
            %% Call position mapping function
            val1 = posmapref(Grids{r}.XPos(i),r+1,'x','+');
            val2 = posmapref(Grids{r}.XPos(i),r+1,'x','-');
            
            %% Assign to Position vectors
                 Grids{r+1}.XPos(1 + count) = val2 ;
                 Grids{r+1}.XPos(1 + count+1) = val1 ;            
            %% Increment counter
            count=count+2;
            
        end
        
        count = 0;
        first_idx = Grids{r}.YInd(3);
        last_idx = Grids{r}.YInd( length(Grids{r}.YInd)-2 );
        for ( j = first_idx : last_idx)
            
            %% Call position mapping function
            val1 = posmapref(Grids{r}.YPos(j),r+1,'y','+');
            val2 = posmapref(Grids{r}.YPos(j),r+1,'y','-');
            
            %% Assign to Position vectors
                 Grids{r+1}.YPos(1 + count) = val2 ;
                 Grids{r+1}.YPos(1 + count+1) = val1 ;
                 
            %% Increment counter
            count=count+2;
            
        end
        
        if (dims == 3)
            count = 0;
            first_idx = Grids{r}.ZInd(3);
            last_idx = Grids{r}.ZInd( length(Grids{r}.ZInd)-2 );
            for ( k = first_idx : last_idx)
                
                %% Call position mapping function
                val1 = posmapref(Grids{r}.ZPos(k),r+1,'z','+');
                val2 = posmapref(Grids{r}.ZPos(k),r+1,'z','-');
                
                %% Assign to Position vectors
                 Grids{r+1}.ZPos(1 + count) = val2 ;
                 Grids{r+1}.ZPos(1 + count+1) = val1 ;
                
                %% Increment counter
                count=count+2;
                
            end
        else
            Grids{2}.ZPos=1; %% 2D default
        end
    end
    
    
    
    %% Assign MACROSCOPIC quantities
    for ( r = 1 : Nref)
        
        %% Get grid sizes
        N_lim = length(Grids{r+1}.XPos);
        M_lim = length(Grids{r+1}.YPos);
        K_lim = length(Grids{r+1}.ZPos);
        
        %% Resize
        Grids{r+1}.u=zeros( N*M*K*dims ,1 );
        Grids{r+1}.rho=zeros( N*M*K ,1 );

        %% Velocity
        LBM_init_vel(r+1);
        
        %% Density
        LBM_init_rho(r+1);
        
    end
    
    %% Generate POPULATION MATRICES for lower levels
    for ( r = 1 : Nref )
        
        %% Get grid sizes
        N_lim = length(Grids{r+1}.XPos);
        M_lim = length(Grids{r+1}.YPos);
        K_lim = length(Grids{r+1}.ZPos);
        
        %% Resize
       	Grids{r+1}.f=zeros(  N_lim * M_lim * K_lim * nVels ,1 );
	    Grids{r+1}.feq=zeros( N_lim * M_lim * K_lim * nVels ,1);
        
        for (i = 1 : length(Grids{r+1}.XInd))
            for (j = 1 : length(Grids{r+1}.YInd))
                for (k = 1 : length(Grids{r+1}.ZInd))
                    for ( v = 1: nVels)
                        
                        
                        %% Initialise f to feq
                        idx  = idxmap(i,j,k,v,M_lim,K_lim,nVels);
                        Grids{r+1}.f(idx) = LBM_collide(i,j,k,v,r+1);
                        
                    end
                end
            end
        end
        Grids{r+1}.feq = Grids{r+1}.f; %% Set feq to feq
        
        %% Compute relaxation time
        Grids{r+1}.omega = 1 / ( ( (1/Grids{r}.omega - .5) *2) + .5);
    end
    
end


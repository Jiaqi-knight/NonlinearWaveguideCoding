% LBM algorithm applicable for both single and multi-grid
% Supply an integer r indicating from which level the algorithm is to be executed.
function LBM_multi (r) 

% Global variables accessible through header file references
global Nref Grids
% Loop twice as refinement ratio per level is 2
count = 1;
while (count < 3)
    
    % Collision on Lr
    if (r == 1)
        
        % Collide on whole grid
        LBM_collide_loop(r, false);
        
    else
        
        % Collide on core only (excludes upper transition layer)
        LBM_collide_loop(r, true);
        
    end
    
    % Check to see if higher level exists and on first loop
    if (r >= 2 && count == 1) 
        
        % Explode and update the upper TL on level r+1
         LBM_explode(r);
        
    end
    
    % Check if lower level exists
    if (Nref+1 > r)
        % Call same routine for lower level
        LBM_multi(r+1);
        
        % Stream
        LBM_stream(r);
        
        
        % Coalesce
        LBM_coalesce(r);
        
    else
        
        % Stream
        LBM_stream(r); %bug
        
        
    end
    
    % Apply boundary conditions
    LBM_boundary(r,0);
    
    % Update macroscopic quantities
    LBM_macro(r);
    
    
    % Check if on L0 and if so drop out as only need to loop once on coarsest level
    if (r == 1)
        break;
    end
    
    % Increment counter
    count=count+1;
    
end


end

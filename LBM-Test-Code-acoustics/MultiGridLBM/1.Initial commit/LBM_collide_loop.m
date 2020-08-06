%Collision operator
%Operates on level r and excludes the upper TL sites if core_flag set to true
function LBM_collide_loop( r,  core_flag )

% 	/*
% 	Loop through the lattice pos to compute the new distribution functions.
% 	Equilibrium based on:
% 	       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4
% 	*/
global Grids nVels dims
%Declarations
N_lim = length(Grids{r}.XPos);
M_lim = length(Grids{r}.YPos);
K_lim = length(Grids{r}.ZPos);
% 	 i_low, j_low, k_low, i_high, j_high, k_high;

%Respond to core flag
if (core_flag)
    %Ignore TL sites to upper level
    i_low = 3; i_high = N_lim-1;
    j_low = 3; j_high = M_lim-1;
    if (dims == 3)
        k_low = 3; k_high = K_lim-1;
    else
        k_low = 1; k_high = K_lim; %if 2D set to default
    end
    
else
    i_low = 1; i_high = N_lim;
    j_low = 1; j_high = M_lim;
    k_low = 1; k_high = K_lim;
end

%Create temporary lattice to prevent overwriting useful populations and initialise with same values as
%pre-collision f grid.
% vector<double> f_new;
% f_new=zeros(length( Grids{r}.f ),1);
%Initialise with current f values
f_new = Grids{r}.f;


%Loop over lattice sites
for ( i = i_low : i_high ) 
    for ( j = j_low : j_high ) 
        for ( k = k_low : k_high ) 
            
            %Get index
            idxijk = idxmap(i,j,k,M_lim,K_lim);
            
            
            %Ignore refined sites
            if (Grids{r}.LatTyp(idxijk) == 2)
                %Do nothing as taken care of on lower level grid
                
            else
                
                
                for ( v = 1 :  nVels)
                    
                    %Get index
                    idxf = idxmap(i,j,k,v,M_lim,K_lim,nVels);
                    
                    %Get feq value by calling overload of collision function
                    Grids{r}.feq(idxf) = LBM_collide(i,j,k,v,r);
                    
                    %Recompute distribution function f
                    f_new(idxf) = (-Grids{r}.omega * (Grids{r}.f(idxf) - Grids{r}.feq(idxf)) ) + Grids{r}.f(idxf);
                    
                end
                
            end
            
        end
    end
end

%Update f from fnew
Grids{r}.f = f_new;


end
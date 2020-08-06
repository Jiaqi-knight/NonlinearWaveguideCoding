%% LatBo.cpp : Defines the entry point for the console application.
clc;clear;close all
%%*	
%%	**************************************************************************************************************
%%	**************************************************************************************************************
%%	**																											**
%%	**												LatBo MAIN													**
%%	**																											**
%%	**************************************************************************************************************
%%	**************************************************************************************************************
%%*
global T deltat dims N M K nVels c w cs Nref rho_in u_0x u_0y u_0z a_x a_y a_z b_x b_y b_z kn km kk

% #include 'stdafx.h'
LBM_definitions
init_globalvars
% #include 'LBM_globalvars.h'			%% Global variable references
% ops_mapping

% using namespace std;	%% Use the standard namespace

%% Entry point

% 	/*
% 	***************************************************************************************************************
% 	********************************************* GENERAL INITIALISE **********************************************
% 	***************************************************************************************************************
% 	*/

	%% Fix output format
	format short; %cout.precision(4);
    global Grids
	%% Timing variables
	% clock_t t_start, t_end; %% Wall clock variables
	t=0;			%% Time step counter, initially zero
	tval = 0;		%% Actual value of physical time (for multi-grid do not necessarily have unit time step)
	totalloops = T/deltat;		%% Total number of loops to be performed (computed)

	

% 	/* ***************************************************************************************************************
% 	*********************************************** LEVEL 0 INITIALISE ***********************************************
% 	*************************************************************************************************************** */
% 	
	%% Time step
	Grids{1}.dt = deltat;

	%% Store spacing
	Lx = b_x - a_x;
	Ly = b_y - a_y;
	Lz = b_z - a_z;
	Grids{1}.dx = 2*(Lx/(2*N));
	Grids{1}.dy = 2*(Ly/(2*M));
	Grids{1}.dz = 2*(Lz/(2*K));

if (dims == 3)
	%% Check that lattice volumes are cubes in 3D
	if ( (Lx/N) ~= (Ly/M) || (Lx/N) ~= (Lz/K) ) 
		disp('Need to have lattice volumes which are cubes -- either change N/M/K or change domain dimensions' );
% 		exit(EXIT_FAILURE);
    end
else 
	%% 2D so need square lattice cells
	if ( (Lx/N) ~= (Ly/M) ) 
		disp ('Need to have lattice cells which are squares -- either change N/M or change domain dimensions' );
% 		exit(EXIT_FAILURE);
    end
end

	%% Refined indices on L0
	Grids{1}.XInd =  Ref_startX+1: Ref_endX+1 ;
	Grids{1}.YInd =  Ref_startY+1: Ref_endY+1 ;
	Grids{1}.ZInd =  Ref_startZ+1: Ref_endZ+1 ;

	%% L0 lattice site coordinates
	Grids{1}.XPos = linspace( a_x + Grids{1}.dx/2, b_x - Grids{1}.dx/2, N );
	Grids{1}.YPos = linspace( a_y + Grids{1}.dy/2, b_y - Grids{1}.dy/2, M );
	Grids{1}.ZPos = linspace( a_z + Grids{1}.dz/2, b_z - Grids{1}.dz/2, K );


	%% Initialise L0 macroscopic quantities
	%% Velocity field
	Grids{1}.u=zeros( N*M*K*dims ,1 );
      LBM_init_vel(1);
%     u=reshape(Grids{1}.u,[N,M,K,dims]);
%     figure;grid off
%     surf(u(:,:,1,2),'LineStyle','none')
%     view(2);
	%% Density field
	Grids{1}.rho=zeros( N*M*K ,1 );
	LBM_init_rho(1);

	%% Initialise L0 matrices (f, feq) and typing matrix
	Grids{1}.f=zeros(  N*M*K*nVels ,1 );
	Grids{1}.feq=zeros( N*M*K*nVels ,1);
	Grids{1}.LatTyp=zeros( N*M*K ,1);

	%% Typing defined as follows:
% 	/*
% 	0 == boundary site
% 	1 == coarse site
% 	2 == fine/refined site
% 	3 == TL to upper (coarser) level
% 	4 == TL to lower (finer) level
% 	*/

	for ( i = 1:N) 
		for ( j = 1:M) 
			for ( k = 1:K) 
				for ( v = 1: nVels) 

					%% Initialise f to feq
					idx  = idxmap(i,j,k,v,M,K,nVels);
					Grids{1}.f(idx) = LBM_collide(i,j,k,v,1);

                end

				%% Label as coarse site
				idx = idxmap(i,j,k,M,K);
				Grids{1}.LatTyp(idx) = 1;

            end
        end
    end
	Grids{1}.feq = Grids{1}.f; %% Make feq = feq too
		
	
	%% Relaxation frequency on L0
	%% Assign relaxation frequency corrected for grid and time step size
	Grids{1}.omega = 1 / ( (nu / (Grids{1}.dt*power(cs,2)) ) + .5 );
	disp (['L0 relaxation time = ',num2str(1/Grids{1}.omega)] );

% 	/* ***************************************************************************************************************
% 	**************************************** REFINED LEVELS INITIALISE ***********************************************
% 	*************************************************************************************************************** */

	if (Nref ~= 0) 

		LBM_init_multi();%ok
		
    end

	disp ('Initialisation Complete...' );

	disp ('Initialising LBM time-stepping...' );

	%% Write out
	disp ('Writing Output to <Grids.out>...' );
% 	for ( r = 1:Nref+1) 
% 		lbm_write3(r,t);
%     end


figure
% 	/* ***************************************************************************************************************
% 	********************************************** LBM PROCEDURE *****************************************************
% 	*************************************************************************************************************** */
% 	
	%% LBM
while (tval < T)
		disp (['\n/%%%%%%%%%%%%%%%% Time Step ', num2str(t+1), ' %%%%%%%%%%%%%%%%/']);

		%% Start the clock
        tic
		%% Call LBM procedure from coarsest level (r = 0)
		LBM_multi(1);%bug

		%% Print Time of loop
        toc

		%% Increment counters
		t=t+1;
		tval =tval+ Grids{1}.dt;
        surf(reshape(Grids{1, 1}.rho,[70,60])  );
        view(2);
        caxis([0.8 1.2])

        pause(0.001)
        
       
end


% 	/* ***************************************************************************************************************
% 	*********************************************** POST PROCESS *****************************************************
% 	*************************************************************************************************************** */

	%% None added yet





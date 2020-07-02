function[]=save3Dtovtk(date,i,itmax,baseoutfilename,latticefolder,meshfolder,flags,N,Nx,Ny,Nz,lattice,vel,m,Fel,Fdamp,Fext,Ekinpdof,Ekinp,Ekin,Epotp,Epot,Dfuncp,Dfunc,covariantbase,contravariantbase,metriccoefficients,g,sqrtg,reciprocalmetriccoefficients,firstChristoffelsymbol,secondChristoffelsymbol,Riemanntensor,Riccitensor,R,rho,u,p,S,gammadot,sigma,mu,tau,cs,Q,feq,f,defgrad,rightCauGretens,EulLagstrain,Cauchystress)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: July 31st, 2014
%
%    Description: 
%          Input: 
%         Output: 
%
%                 Ordering of flags:
%                 flags(1)  --> save data over fixed lattice (computational domain)
%                 flags(2)  --> save data over moving mesh (physical domain)
%                 flags(3)  --> save velocities (LSM)
%                 flags(4)  --> save masses (LSM)
%                 flags(5)  --> save elastic forces (LSM)
%                 flags(6)  --> save damping forces (LSM)
%                 flags(7)  --> save external forces (LSM)
%                 flags(8)  --> save total forces (LSM)
%                 flags(9)  --> save particle kinetic energy per dof (LSM)
%                 flags(10) --> save particle kinetic energy (LSM)
%                 flags(11) --> save total kinetic energy (LSM)
%                 flags(12) --> save particle potential energy (LSM)
%                 flags(13) --> save total potential energy (LSM)
%                 flags(14) --> save particle dissipation function (LSM)
%                 flags(15) --> save total dissipation function (LSM)
%                 flags(16) --> save covariant base vectors
%                 flags(17) --> save contravariant base vectors
%                 flags(18) --> save metric tensor
%                 flags(19) --> save metric tensor determinant
%                 flags(20) --> save g
%                 flags(21) --> save reciprocal metric coefficients
%                 flags(22) --> save first Christoffel symbols
%                 flags(23) --> save second Christoffel symbols
%                 flags(24) --> save Riemann tensor
%                 flags(25) --> save Ricci tensor
%                 flags(26) --> save Ricci curvature
%                 flags(27) --> save density (LBM)
%                 flags(28) --> save velocity in cartesian components (LBM)
%                 flags(29) --> save pressure (LBM)
%                 flags(30) --> save strain rate tensor (LBM)
%                 flags(31) --> save strain rate (LBM)
%                 flags(32) --> save stress tensor (LBM)
%                 flags(33) --> save viscosity (LBM)
%                 flags(34) --> save relaxation time (LBM)
%                 flags(35) --> save sound velocity (LBM)
%                 flags(36) --> save particle equilibrium populations (LBM)
%                 flags(37) --> save particle populations (LBM)
%                 flags(38) --> save deformation gradient (LSM)
%                 flags(39) --> save right Cauchy-Green tensor (LSM)
%                 flags(40) --> save Euler-Lagrange strain tensor (LSM)
%                 flags(41) --> save Cauchy stress tensor (LSM)

%%

if i~=0
    idigits = fix(abs(log10(abs(i))))+1;
else
    idigits = 1;
end
if itmax~=0
    itmaxdigits = fix(abs(log10(abs(itmax))))+1;
else
    itmaxdigits = 1;
end
if idigits<itmaxdigits
    itstring = num2str(i);
    diffdigits = itmaxdigits - idigits;
    for i=1:diffdigits
        itstring = strcat('0',itstring);
    end
else
    itstring = '00';
end

Ncells = (Nx-1)*(Ny-1)*(Nz-1);
cells = generatecells3D(Nx,Ny,Nz);
cells(:,2:9) = cells(:,2:9) - 1;
   
if flags(1)
    % Create file

    latticefilename = strcat(latticefolder,'/',baseoutfilename,'_N',itstring,'.vtk');

    latticefid = fopen(latticefilename,'w');

    % Header
    % lines: 5
    % The first line simply describes the version number. The second line can contain what ever you want, for example a bit of information about who generated the file. The third simply determines the filetype: ASCII. The forth describes what kind of data VTK can expected in the rest of the file. There are five different types: Three structured types (STRUCTURED_POINTS, STRUCTURED_GRID, RECTILINEAR_GRID) in which all points somehow lineup or match, not very usable. One surface type (POLYDATA) which doesn't contain volume, not what we need. And finally UNSTRUCTURED_GRID which describes arbitrary combinations of all availible cells types. This is what we like.
    % The fifth line is empty.
    % --> Unstructured grid: 2D or 3D grid; for every grid point all three coordinates and for each grid cell all constituent points and the cell shape are given. XML file extension: *.vtu
    % --> Polygonal data:    2D grid; like the unstructured grid, but there are no polyhedra, but only flat polygons. Especially suited for maps (topography). XML file extension: *.vtp
    % --> Rectilinear grid:  3D grid; the grid cells are cuboids, so only the steps along the coordinate axes have to be given, but not the individual point coordinates or the connectivity.
    % --> Structured grid:   3D grid; here, all point coordinates are given, but the connectivity is omitted.
    % --> Structured points: Like the rectilinear grid, but the spacing between the points is equidistant; so only the origin and the spacing has to be given, not the point coordinates.

    fprintf(latticefid,'# vtk DataFile Version 6.0\n');
    fprintf(latticefid,strcat('Data over lattice in computational domain. Generated by Luca Di Stasio, Computational Physics of Engineering Materials, Institute for Building Materials, ETH Zuerich on ',date,'\n'));
    fprintf(latticefid,'ASCII\n');
    fprintf(latticefid,'DATASET STRUCTURED_GRID');
    fprintf(latticefid,' \n');

    % Pointdata
    % lines: #nr_points+2
    % The first line notes how many points there will be ("9") and in which format they'll be supplied (usually "FLOAT"). Nine lines follow, each line containing the xyz-coordinates of a point.
    % Based on this list, each point is assigned an id. The first point has id 0, the second point id 1, and so forth, till the last point with id 8.
    % The last line is empty.

    fprintf(latticefid,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
    fprintf(latticefid,'POINTS %d FLOAT\n',N);
    for i=1:N
        fprintf(latticefid,'%12.8f %12.8f %12.8f\n',lattice(i,1),lattice(i,2),lattice(i,3));
    end
    fprintf(latticefid,' \n');
    
    % Celldata
    % lines: #nr_cells+5
    % The first line notes how many cells there will be ("4") and how many numbers total will be supplied in the CELLS-block ("20"). Four lines follow, each line containing the information for a cell.
    % Each cell line starts with a number saying how many point-ids are to be read in that line followed by the list of those point-ids.
    % Based on this list, each cell is assigned an id. The first cell has id 0, the second id 1, and so forth, till the last cell with id 3.
    % After all cells have been given points, the cell types have to be set.
    % This is done in the CELL_TYPES-block.
    % The CELL_TYPES-line says how many cell-types are to be set. This number has to be the same as the number in the CELLS-line ("4"). The following line is a list of cell-types that are assigned to all cells. The first cell is of type "9", so ist the second cell and so forth.
    % The last line is empty.

    % Pointdatasets
    % lines: #nr_datasets*(#nr_points +3)+1
    % In this segment of the file, datasets that are associated to the points of the grid can be given. The data can either be of scalar- or vector-type.
    % The first line initiates that pointdata will be supplied in the following lines, the "9" simply says how many points per dataset there are (yes, it's redundant).
    % We seperate each of the datasets with an empty line. Each of the datasets has two leading lines followed by "9" lines of data.
    % The first dataset line specifies if "SCALAR"- or "VECTOR"-data will follow. It also assigns a name to the dataset ("HorizontalSpeed" , "Temperature"). The "FLOAT" specifies the format of the data. "FLOAT" will most likely always be good enough.
    % The second line has to do with color tables and "LOOKUP_TABLE default" will most likely be the only setting you will ever need.
    % The following "9" lines will either contain per line a scalar-value or a vector-value (consisting of 3 numbers).
    % The last line is empty.

    if any(flags(3:end)==1)
        fprintf(latticefid,'POINT_DATA %d\n',N);
    end
    if flags(3)
        fprintf(latticefid,'VECTORS Velocity FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',vel(i,1),vel(i,2),vel(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(4)
        fprintf(latticefid,'SCALARS Particles_mass FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',m(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(5)
        fprintf(latticefid,'VECTORS Elastic_forces FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Fel(i,1),Fel(i,2),Fel(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(6)
        fprintf(latticefid,'VECTORS Damping_forces FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Fdamp(i,1),Fdamp(i,2),Fdamp(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(7)
        fprintf(latticefid,'VECTORS External_forces FLOAT\n');
        Fext = full(Fext);
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Fext(i,1),Fext(i,2),Fext(i,3));
        end
        fprintf(latticefid,' \n');
        Fext = sparse(Fext);
    end
    if flags(8)
        fprintf(latticefid,'VECTORS Total_forces FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Fel(i,1)+Fdamp(i,1)+Fext(i,1),Fel(i,2)+Fdamp(i,2)+Fext(i,2),Fel(i,3)+Fdamp(i,3)+Fext(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(9)
        fprintf(latticefid,'VECTORS Particles_kinetic_energy_dof FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Ekinpdof(i,1),Ekinpdof(i,2),Ekinpdof(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(10)
        fprintf(latticefid,'SCALARS Particles_kinetic_energy FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Ekinp(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(11)
        fprintf(latticefid,'SCALARS Kinetic_energy FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Ekin);
        end
        fprintf(latticefid,' \n');
    end
    if flags(12)
        fprintf(latticefid,'SCALARS Particles_potential_energy FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Epotp(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(13)
        fprintf(latticefid,'SCALARS Potential_energy FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Epot);
        end
        fprintf(latticefid,' \n');
    end
    if flags(14)
        fprintf(latticefid,'SCALARS Particles_dissipation_function FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Dfuncp(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(15)
        fprintf(latticefid,'SCALARS Dissipation_function FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',Dfunc);
        end
        fprintf(latticefid,' \n');
    end
    if flags(16)
        fprintf(latticefid,'VECTORS Covariant_base FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',covariantbase(i,1),covariantbase(i,2),covariantbase(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(17)
        fprintf(latticefid,'VECTORS Contravariant_base FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',contravariantbase(i,1),contravariantbase(i,2),contravariantbase(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(18)
        fprintf(latticefid,'TENSORS Metric_coefficients FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,1),metriccoefficients(i,4),metriccoefficients(i,5));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,4),metriccoefficients(i,2),metriccoefficients(i,6));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,5),metriccoefficients(i,6),metriccoefficients(i,3));
            fprintf(latticefid,' \n');
        end
    end
    if flags(19)
        fprintf(latticefid,'SCALARS Metric_determinant_J FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',g(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(20)
        fprintf(latticefid,'SCALARS Sqrt_metric_determinant_g FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',sqrtg(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(21)
        fprintf(latticefid,'TENSORS Reciprocal_metric_coefficients FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,1),reciprocalmetriccoefficients(i,4),reciprocalmetriccoefficients(i,5));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,4),reciprocalmetriccoefficients(i,2),reciprocalmetriccoefficients(i,6));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,5),reciprocalmetriccoefficients(i,6),reciprocalmetriccoefficients(i,3));
            fprintf(latticefid,' \n');
        end
    end
    if flags(22)
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_111 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,1));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_121 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,2));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_131 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,3));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_211 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,4));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_221 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,5));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_231 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,6));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_311 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,7));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_321 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,8));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_331 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,9));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_112 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,10));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_122 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,11));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_132 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,12));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_212 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,13));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_222 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,14));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_232 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,15));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_312 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,16));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_322 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,17));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_332 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,18));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_113 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,19));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_123 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,20));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_133 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,21));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_213 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,22));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_223 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,23));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_233 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,24));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_313 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,25));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_323 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,26));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS First_Christoffel_symbol_component_333 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',firstChristoffelsymbol(i,27));
        end
        fprintf(latticefid,' \n');
    end
    if flags(23)
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_111 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,1));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_121 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,2));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_131 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,3));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_211 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,4));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_221 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,5));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_231 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,6));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_311 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,7));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_321 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,8));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_331 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,9));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_112 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,10));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_122 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,11));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_132 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,12));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_212 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,13));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_222 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,14));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_232 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,15));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_312 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,16));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_322 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,17));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_332 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,18));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_113 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,19));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_123 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,20));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_133 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,21));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_213 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,22));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_223 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,23));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_233 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,24));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_313 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,25));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_323 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,26));
        end
        fprintf(latticefid,' \n');
        fprintf(latticefid,'SCALARS Second_Christoffel_symbol_component_333 FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',secondChristoffelsymbol(i,27));
        end
        fprintf(latticefid,' \n');
    end
    if flags(24)
        Riemanntensorindices = [1111 1112 1113 1121 1122 1123 1131 1132 1133 1211 1212 1213 1221 1222 1223 1231 1232 1233 1311 1312 1313 1321 1322 1323 1331 1332 1333 2111 2112 2113 2121 2122 2123 2131 2132 2133 2211 2212 2213 2221 2222 2223 2231 2232 2233 2311 2312 2313 2321 2322 2323 2331 2332 2333 3111 3112 3113 3121 3122 3123 3131 3132 3133 3211 3212 3213 3221 3222 3223 3231 3232 3233 3311 3312 3313 3321 3322 3323 3331 3332 3333];
        for j=1:length(Riemanntensorindices)
            fprintf(latticefid,strcat('SCALARS Riemann_tensor_component_',num2str(Riemanntensorindices(j)),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',Riemanntensor(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(25)
        fprintf(latticefid,'TENSORS Ricci_tensor FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,1),Riccitensor(i,4),Riccitensor(i,5));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,4),Riccitensor(i,2),Riccitensor(i,6));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,5),Riccitensor(i,6),Riccitensor(i,3));
            fprintf(latticefid,' \n');
        end
    end
    if flags(26)
        fprintf(latticefid,'SCALARS Ricci_scalar FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',R(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(27)
        fprintf(latticefid,'SCALARS Fluid_density FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',rho(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(28)
        fprintf(latticefid,'VECTORS Fluid_velocity FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',u(i,1),u(i,2),u(i,3));
        end
        fprintf(latticefid,' \n');
    end
    if flags(29)
        fprintf(latticefid,'SCALARS Fluid_pressure FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',p(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(30)
        fprintf(latticefid,'TENSORS Fluid_strain_rate_tensor FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',S(i,1),S(i,4),S(i,5));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',S(i,4),S(i,2),S(i,6));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',S(i,5),S(i,6),S(i,3));
            fprintf(latticefid,' \n');
        end
    end
    if flags(31)
        fprintf(latticefid,'SCALARS Fluid_strain_rate FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',gammadot(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(32)
        fprintf(latticefid,'TENSORS Fluid_stress_tensor FLOAT\n');
        for i=1:N
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',sigma(i,1),sigma(i,4),sigma(i,5));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',sigma(i,4),sigma(i,2),sigma(i,6));
            fprintf(latticefid,'%12.8f %12.8f %12.8f\n',sigma(i,5),sigma(i,6),sigma(i,3));
            fprintf(latticefid,' \n');
        end
    end
    if flags(33)
        fprintf(latticefid,'SCALARS Fluid_viscosity FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',mu(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(34)
        fprintf(latticefid,'SCALARS Fluid_relaxation_time FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',tau(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(35)
        fprintf(latticefid,'SCALARS Fluid_sound_velocity FLOAT\n');
        fprintf(latticefid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(latticefid,'%12.8f\n',cs(i));
        end
        fprintf(latticefid,' \n');
    end
    if flags(36)
        for j=1:Q
            fprintf(latticefid,strcat('SCALARS Fluid_particle_equilibrium_populations_Q',num2str(j),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',feq(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(37)
        for j=1:Q
            fprintf(latticefid,strcat('SCALARS Fluid_particle_populations_Q',num2str(j),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',f(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(38)
        defgradindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(defgradindices)
            fprintf(latticefid,strcat('SCALARS Deformation_gradient_component_',num2str(defgradindices(j)),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',defgrad(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(39)
        rightCauGretensindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(rightCauGretensindices)
            fprintf(latticefid,strcat('SCALARS Right_Cauchy_Green_tensor_component_',num2str(rightCauGretensindices(j)),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',rightCauGretens(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(40)
        EulLagstrainindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(EulLagstrainindices)
            fprintf(latticefid,strcat('SCALARS Euler_Lagrange_strain_tensor_component_',num2str(EulLagstrainindices(j)),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',EulLagstrain(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end
    if flags(41)
        Cauchystressindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(Cauchystressindices)
            fprintf(latticefid,strcat('SCALARS Cauchy_stress_tensor_component_',num2str(Cauchystressindices(j)),' FLOAT\n'));
            fprintf(latticefid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(latticefid,'%12.8f\n',Cauchystress(i,j));
            end
            fprintf(latticefid,' \n');
        end
    end

    % Celldatasets
    % lines: #nr_datasets*(#nr_cells +3)+1
    % Exactly the same setup as the point-data, only this data is assigned to cells. 

    % Close file

    fclose(latticefid);
end

if flags(2)
    % Create file
    
    meshfilename = strcat(meshfolder,'/',baseoutfilename,'_N',itstring,'.vtk');

    meshfid = fopen(meshfilename,'w');

    % Header
    % lines: 5
    % The first line simply describes the version number. The second line can contain what ever you want, for example a bit of information about who generated the file. The third simply determines the filetype: ASCII. The forth describes what kind of data VTK can expected in the rest of the file. There are five different types: Three structured types (STRUCTURED_POINTS, STRUCTURED_GRID, RECTILINEAR_GRID) in which all points somehow lineup or match, not very usable. One surface type (POLYDATA) which doesn't contain volume, not what we need. And finally UNSTRUCTURED_GRID which describes arbitrary combinations of all availible cells types. This is what we like.
    % The fifth line is empty.
    % --> Unstructured grid: 2D or 3D grid; for every grid point all three coordinates and for each grid cell all constituent points and the cell shape are given. XML file extension: *.vtu
    % --> Polygonal data:    2D grid; like the unstructured grid, but there are no polyhedra, but only flat polygons. Especially suited for maps (topography). XML file extension: *.vtp
    % --> Rectilinear grid:  3D grid; the grid cells are cuboids, so only the steps along the coordinate axes have to be given, but not the individual point coordinates or the connectivity.
    % --> Structured grid:   3D grid; here, all point coordinates are given, but the connectivity is omitted.
    % --> Structured points: Like the rectilinear grid, but the spacing between the points is equidistant; so only the origin and the spacing has to be given, not the point coordinates.

    fprintf(meshfid,'# vtk DataFile Version 6.0\n');
    fprintf(meshfid,strcat('Data over mesh in physical domain. Generated by Luca Di Stasio, Computational Physics of Engineering Materials, Institute for Building Materials, ETH Zuerich on ',date,'\n'));
    fprintf(meshfid,'ASCII\n');
    fprintf(meshfid,'DATASET UNSTRUCTURED_GRID');
    fprintf(meshfid,' \n');

    % Pointdata
    % lines: #nr_points+2
    % The first line notes how many points there will be ("9") and in which format they'll be supplied (usually "FLOAT"). Nine lines follow, each line containing the xyz-coordinates of a point.
    % Based on this list, each point is assigned an id. The first point has id 0, the second point id 1, and so forth, till the last point with id 8.
    % The last line is empty.

    fprintf(meshfid,'POINTS %d FLOAT\n',N);
    for i=1:N
        fprintf(meshfid,'%12.8f %12.8f %12.8f\n',lattice(i,7),lattice(i,8),lattice(i,9));
    end
    fprintf(meshfid,' \n');

    % Celldata
    % lines: #nr_cells+5
    % The first line notes how many cells there will be ("4") and how many numbers total will be supplied in the CELLS-block ("20"). Four lines follow, each line containing the information for a cell.
    % Each cell line starts with a number saying how many point-ids are to be read in that line followed by the list of those point-ids.
    % Based on this list, each cell is assigned an id. The first cell has id 0, the second id 1, and so forth, till the last cell with id 3.
    % After all cells have been given points, the cell types have to be set.
    % This is done in the CELL_TYPES-block.
    % The CELL_TYPES-line says how many cell-types are to be set. This number has to be the same as the number in the CELLS-line ("4"). The following line is a list of cell-types that are assigned to all cells. The first cell is of type "9", so ist the second cell and so forth.
    % The last line is empty.

    fprintf(meshfid,'CELLS %d %d\n',Ncells,Ncells*9);
    for i=1:Ncells
        fprintf(meshfid,'%d %d %d %d %d %d %d %d %d\n',cells(i,1),cells(i,2),cells(i,3),cells(i,4),cells(i,5),cells(i,6),cells(i,7),cells(i,8),cells(i,9));
    end
    fprintf(meshfid,' \n');

    fprintf(meshfid,'CELL_TYPES %d\n',Ncells);
    for i=1:Ncells
        fprintf(meshfid,'%d\n',12);
    end
    fprintf(meshfid,' \n');

    % Pointdatasets
    % lines: #nr_datasets*(#nr_points +3)+1
    % In this segment of the file, datasets that are associated to the points of the grid can be given. The data can either be of scalar- or vector-type.
    % The first line initiates that pointdata will be supplied in the following lines, the "9" simply says how many points per dataset there are (yes, it's redundant).
    % We seperate each of the datasets with an empty line. Each of the datasets has two leading lines followed by "9" lines of data.
    % The first dataset line specifies if "SCALAR"- or "VECTOR"-data will follow. It also assigns a name to the dataset ("HorizontalSpeed" , "Temperature"). The "FLOAT" specifies the format of the data. "FLOAT" will most likely always be good enough.
    % The second line has to do with color tables and "LOOKUP_TABLE default" will most likely be the only setting you will ever need.
    % The following "9" lines will either contain per line a scalar-value or a vector-value (consisting of 3 numbers).
    % The last line is empty.

    if any(flags(3:end)==1)
        fprintf(meshfid,'POINT_DATA %d\n',N);
    end
    if flags(3)
        fprintf(meshfid,'VECTORS Velocity FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',vel(i,1),vel(i,2),vel(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(4)
        fprintf(meshfid,'SCALARS Particles_mass FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',m(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(5)
        fprintf(meshfid,'VECTORS Elastic_forces FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Fel(i,1),Fel(i,2),Fel(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(6)
        fprintf(meshfid,'VECTORS Damping_forces FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Fdamp(i,1),Fdamp(i,2),Fdamp(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(7)
        fprintf(meshfid,'VECTORS External_forces FLOAT\n');
        Fext = full(Fext);
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Fext(i,1),Fext(i,2),Fext(i,3));
        end
        fprintf(meshfid,' \n');
        Fext = sparse(Fext);
    end
    if flags(8)
        fprintf(meshfid,'VECTORS Total_forces FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Fel(i,1)+Fdamp(i,1)+Fext(i,1),Fel(i,2)+Fdamp(i,2)+Fext(i,2),Fel(i,3)+Fdamp(i,3)+Fext(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(9)
        fprintf(meshfid,'VECTORS Particles_kinetic_energy_dof FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Ekinpdof(i,1),Ekinpdof(i,2),Ekinpdof(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(10)
        fprintf(meshfid,'SCALARS Particles_kinetic_energy FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Ekinp(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(11)
        fprintf(meshfid,'SCALARS Kinetic_energy FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Ekin);
        end
        fprintf(meshfid,' \n');
    end
    if flags(12)
        fprintf(meshfid,'SCALARS Particles_potential_energy FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Epotp(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(13)
        fprintf(meshfid,'SCALARS Potential_energy FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Epot);
        end
        fprintf(meshfid,' \n');
    end
    if flags(14)
        fprintf(meshfid,'SCALARS Particles_dissipation_function FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Dfuncp(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(15)
        fprintf(meshfid,'SCALARS Dissipation_function FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',Dfunc);
        end
        fprintf(meshfid,' \n');
    end
    if flags(16)
        fprintf(meshfid,'VECTORS Covariant_base FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',covariantbase(i,1),covariantbase(i,2),covariantbase(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(17)
        fprintf(meshfid,'VECTORS Contravariant_base FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',contravariantbase(i,1),contravariantbase(i,2),contravariantbase(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(18)
        fprintf(meshfid,'TENSORS Metric_coefficients FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,1),metriccoefficients(i,4),metriccoefficients(i,5));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,4),metriccoefficients(i,2),metriccoefficients(i,6));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',metriccoefficients(i,5),metriccoefficients(i,6),metriccoefficients(i,3));
            fprintf(meshfid,' \n');
        end
    end
    if flags(19)
        fprintf(meshfid,'SCALARS Metric_determinant_J FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',g(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(20)
        fprintf(meshfid,'SCALARS Sqrt_metric_determinant_g FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',sqrtg(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(21)
        fprintf(meshfid,'TENSORS Reciprocal_metric_coefficients FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,1),reciprocalmetriccoefficients(i,4),reciprocalmetriccoefficients(i,5));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,4),reciprocalmetriccoefficients(i,2),reciprocalmetriccoefficients(i,6));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',reciprocalmetriccoefficients(i,5),reciprocalmetriccoefficients(i,6),reciprocalmetriccoefficients(i,3));
            fprintf(meshfid,' \n');
        end
    end
    if flags(22)
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_111 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,1));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_121 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,2));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_131 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,3));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_211 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,4));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_221 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,5));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_231 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,6));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_311 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,7));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_321 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,8));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_331 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,9));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_112 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,10));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_122 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,11));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_132 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,12));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_212 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,13));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_222 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,14));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_232 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,15));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_312 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,16));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_322 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,17));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_332 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,18));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_113 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,19));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_123 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,20));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_133 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,21));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_213 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,22));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_223 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,23));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_233 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,24));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_313 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,25));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_323 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,26));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS First_Christoffel_symbol_component_333 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',firstChristoffelsymbol(i,27));
        end
        fprintf(meshfid,' \n');
    end
    if flags(23)
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_111 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,1));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_121 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,2));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_131 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,3));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_211 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,4));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_221 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,5));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_231 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,6));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_311 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,7));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_321 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,8));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_331 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,9));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_112 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,10));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_122 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,11));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_132 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,12));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_212 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,13));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_222 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,14));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_232 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,15));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_312 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,16));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_322 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,17));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_332 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,18));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_113 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,19));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_123 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,20));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_133 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,21));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_213 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,22));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_223 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,23));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_233 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,24));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_313 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,25));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_323 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,26));
        end
        fprintf(meshfid,' \n');
        fprintf(meshfid,'SCALARS Second_Christoffel_symbol_component_333 FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',secondChristoffelsymbol(i,27));
        end
        fprintf(meshfid,' \n');
    end
    if flags(24)
        Riemanntensorindices = [1111 1112 1113 1121 1122 1123 1131 1132 1133 1211 1212 1213 1221 1222 1223 1231 1232 1233 1311 1312 1313 1321 1322 1323 1331 1332 1333 2111 2112 2113 2121 2122 2123 2131 2132 2133 2211 2212 2213 2221 2222 2223 2231 2232 2233 2311 2312 2313 2321 2322 2323 2331 2332 2333 3111 3112 3113 3121 3122 3123 3131 3132 3133 3211 3212 3213 3221 3222 3223 3231 3232 3233 3311 3312 3313 3321 3322 3323 3331 3332 3333];
        for j=1:length(Riemanntensorindices)
            fprintf(meshfid,strcat('SCALARS Riemann_tensor_component_',num2str(Riemanntensorindices(j)),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',Riemanntensor(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(25)
        fprintf(meshfid,'TENSORS Ricci_tensor FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,1),Riccitensor(i,4),Riccitensor(i,5));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,4),Riccitensor(i,2),Riccitensor(i,6));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',Riccitensor(i,5),Riccitensor(i,6),Riccitensor(i,3));
            fprintf(meshfid,' \n');
        end
    end
    if flags(26)
        fprintf(meshfid,'SCALARS Ricci_scalar FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',R(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(27)
        fprintf(meshfid,'SCALARS Fluid_density FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',rho(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(28)
        fprintf(meshfid,'VECTORS Fluid_velocity FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',u(i,1),u(i,2),u(i,3));
        end
        fprintf(meshfid,' \n');
    end
    if flags(29)
        fprintf(meshfid,'SCALARS Fluid_pressure FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',p(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(30)
        fprintf(meshfid,'TENSORS Fluid_strain_rate_tensor FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',S(i,1),S(i,4),S(i,5));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',S(i,4),S(i,2),S(i,6));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',S(i,5),S(i,6),S(i,3));
            fprintf(meshfid,' \n');
        end
    end
    if flags(31)
        fprintf(meshfid,'SCALARS Fluid_strain_rate FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',gammadot(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(32)
        fprintf(meshfid,'TENSORS Fluid_stress_tensor FLOAT\n');
        for i=1:N
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',sigma(i,1),sigma(i,4),sigma(i,5));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',sigma(i,4),sigma(i,2),sigma(i,6));
            fprintf(meshfid,'%12.8f %12.8f %12.8f\n',sigma(i,5),sigma(i,6),sigma(i,3));
            fprintf(meshfid,' \n');
        end
    end
    if flags(33)
        fprintf(meshfid,'SCALARS Fluid_viscosity FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',mu(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(34)
        fprintf(meshfid,'SCALARS Fluid_relaxation_time FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',tau(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(35)
        fprintf(meshfid,'SCALARS Fluid_sound_velocity FLOAT\n');
        fprintf(meshfid,'LOOKUP_TABLE default\n');
        for i=1:N
            fprintf(meshfid,'%12.8f\n',cs(i));
        end
        fprintf(meshfid,' \n');
    end
    if flags(36)
        for j=1:Q
            fprintf(meshfid,strcat('SCALARS Fluid_particle_equilibrium_populations_Q',num2str(j),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',feq(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(37)
        for j=1:Q
            fprintf(meshfid,strcat('SCALARS Fluid_particle_populations_Q',num2str(j),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',f(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(38)
        defgradindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(defgradindices)
            fprintf(meshfid,strcat('SCALARS Deformation_gradient_component_',num2str(defgradindices(j)),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',defgrad(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(39)
        rightCauGretensindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(rightCauGretensindices)
            fprintf(meshfid,strcat('SCALARS Right_Cauchy_Green_tensor_component_',num2str(rightCauGretensindices(j)),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',rightCauGretens(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(40)
        EulLagstrainindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(EulLagstrainindices)
            fprintf(meshfid,strcat('SCALARS Euler_Lagrange_strain_tensor_component_',num2str(EulLagstrainindices(j)),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',EulLagstrain(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end
    if flags(41)
        Cauchystressindices = [11 12 13 21 22 23 31 32 33];
        for j=1:length(Cauchystressindices)
            fprintf(meshfid,strcat('SCALARS Cauchy_stress_tensor_component_',num2str(Cauchystressindices(j)),' FLOAT\n'));
            fprintf(meshfid,'LOOKUP_TABLE default\n');
            for i=1:N
                fprintf(meshfid,'%12.8f\n',Cauchystress(i,j));
            end
            fprintf(meshfid,' \n');
        end
    end

    % Celldatasets
    % lines: #nr_datasets*(#nr_cells +3)+1
    % Exactly the same setup as the point-data, only this data is assigned to cells. 

    % Close file

    fclose(meshfid);

end

return

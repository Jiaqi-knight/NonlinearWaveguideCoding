function[lattice]=generate_lattice(D,dim1min,deltadim1,Ndim1,dim2min,deltadim2,Ndim2,dim3min,deltadim3,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 1st, 2014
%    Last update: June 12th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

lattice.D = D;
lattice.dim1min = dim1min;
lattice.dim1max = dim1min + deltadim1*Ndim1;
lattice.Ndim1 = Ndim1;
lattice.dim2min = dim2min;
lattice.dim2max = dim2min + deltadim2*Ndim2;
lattice.Ndim2 = Ndim2;
lattice.dim3min = dim3min;
lattice.dim3max = dim3min + deltadim3*Ndim3;
lattice.Ndim3 = Ndim3;
lattice.Nlinesdim1 = Nlinesdim1;
lattice.Nlinesdim2 = Nlinesdim2;
lattice.Nlinesdim3 = Nlinesdim3;
lattice.totN = Ndim1*Ndim2*Ndim3;
lattice.edgeflag = edgeflag;
lattice.faceflag = faceflag;
lattice.cellflag = cellflag;
   
switch D
%-----------------------------------------------------------------------------1D--------------------------------------------------------------------------------------------
    case 1
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('lattice GENERATION')
            disp(' ')
            disp('Generation of equally spaced intervals along a line - 1D')
            fprintf('Start point x0 = %d\n',lattice.dim1min)
            fprintf('End point xN = %6.2f\n',lattice.dim1max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',lattice.totN)
        end
        lattice.nodes.id = zeros(1,lattice.totN);
        lattice.nodes.meshcoordinates = zeros(3,lattice.totN);  
        lattice.deltadim1 = deltadim1;
        for i=1:Ndim1
            lattice.nodes.meshcoordinates(1,i) = dim1min + deltadim1*(i-1);
            lattice.nodes.meshcoordinates(2,i) = 0; 
            lattice.nodes.meshcoordinates(3,i) = 0;
            lattice.nodes.id(i) = i;
        end
        lattice.coordlines.dim1free = zeros(lattice.Nlinesdim1,3);
        coord2 = 0;
        coord3 = 0;
        deltalinedim1 = (lattice.dim1max - lattice.dim1min)/(lattice.Nlinesdim1 - 1);
        for l=1:lattice.Nlinesdim1
            lattice.coordlines.dim1free(l,1) = lattice.dim1max + deltalinedim1*(l-1);
            lattice.coordlines.dim1free(l,2) = coord2; 
            lattice.coordlines.dim1free(l,3) = coord3;
        end
        if edgeflag || faceflag || cellflag
            lattice.edgeflag = 1;
            lattice.edgespernode = 2;
            lattice.nodesperedge = 2;
            lattice.totE = (lattice.Ndim1-1);
            lattice.nodes.edges = zeros(lattice.edgespernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(lattice.nodesperedge,lattice.totE);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
            end
            for i=1:lattice.Ndim1-1
               nodeindex1 = i;
               nodeindex2 = i + 1;
               edgeindex1 = i;
               lattice.nodes.edges(2,nodeindex1) = edgeindex1;
               lattice.nodes.edges(1,nodeindex2) = edgeindex1;
               lattice.edges.id(i) = edgeindex1;
               lattice.edges.centroid(1,i) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
               lattice.edges.centroid(2,i) = 0;
               lattice.edges.centroid(3,i) = 0;
               lattice.edges.nodes(1,i) = nodeindex1;
               lattice.edges.nodes(2,i) = nodeindex2;
            end
        end
%-----------------------------------------------------------------------------2D--------------------------------------------------------------------------------------------
    case 2
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('lattice GENERATION')
            disp(' ')
            disp('Meshing of a rectangular domain with quadrilateral elements - 2D')
            fprintf('Lower boundary along the first direction x0 = %d\n',lattice.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',lattice.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',lattice.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',lattice.dim2max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',lattice.totN)
        end
        lattice.nodes.id = zeros(1,lattice.totN);
        lattice.nodes.meshcoordinates = zeros(3,lattice.totN);  
        lattice.deltadim1 = deltadim1;
        lattice.deltadim2 = deltadim2;
        for i=1:Ndim2
            for j=1:Ndim1
                lattice.nodes.meshcoordinates(1,j + (i-1)*Ndim1) = dim1min + deltadim1*(j-1);
                lattice.nodes.meshcoordinates(2,j + (i-1)*Ndim1) = dim2min + deltadim2*(i-1); 
                lattice.nodes.meshcoordinates(3,j + (i-1)*Ndim1) = 0;
                lattice.nodes.id(j + i*Ndim1) = j + (i-1)*Ndim1;
            end
        end
        lattice.coordlines.dim1free = zeros(lattice.Nlinesdim1,3*lattice.Ndim2);
        lattice.coordlines.dim2free = zeros(lattice.Nlinesdim2,3*lattice.Ndim1);
        for j=1:lattice.Ndim2
            coord2 = lattice.dim2min + lattice.deltadim2*(j-1);
            coord3 = 0;
            deltalinedim1 = (lattice.dim1max - lattice.dim1min)/(lattice.Nlinesdim1 - 1);
            for l=1:lattice.Nlinesdim1
                lattice.coordlines.dim1free(l,3*(j-1) + 1) = lattice.dim1max + deltalinedim1*(l-1);
                lattice.coordlines.dim1free(l,3*(j-1) + 2) = coord2; 
                lattice.coordlines.dim1free(l,3*(j-1) + 3) = coord3;
            end
        end
        for k=1:lattice.Ndim1
            coord1 = lattice.dim1min + lattice.deltadim1*(k-1);
            coord3 = 0;
            deltalinedim2 = (lattice.dim2max - lattice.dim2min)/(lattice.Nlinesdim2 - 1);
            for l=1:lattice.Nlinesdim2
                lattice.coordlines.dim2free(l,3*(k-1) + 1) = coord1;
                lattice.coordlines.dim2free(l,3*(k-1) + 2) = lattice.dim2max + deltalinedim2*(l-1); 
                lattice.coordlines.dim2free(l,3*(k-1) + 3) = coord3;
            end
        end
        if (edgeflag && faceflag) || (edgeflag && cellflag)
            lattice.edgeflag = 1;
            lattice.faceflag = 1;
            lattice.totE = (lattice.Ndim1-1)*lattice.Ndim2 + (lattice.Ndim2-1)*lattice.Ndim1;
            lattice.totF = (lattice.Ndim1-1)*(lattice.Ndim2-1);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
                fprintf('Total number of faces F = %6.2f\n',lattice.totF)
                disp(' ')
            end
            lattice.nodesperedge = 2;
            lattice.nodesperface = 4;
            lattice.edgespernode = 4;
            lattice.edgesperface = 4;
            lattice.facespernode = 4;
            lattice.facesperedge = 2;
            lattice.nodes.edges = zeros(edgespernode,lattice.totN);
            lattice.nodes.faces = zeros(facespernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(nodesperedge,lattice.totE);
            lattice.edges.faces = zeros(facesperedge,lattice.totE);
            lattice.faces.id = zeros(1,lattice.totF);
            lattice.faces.centroid = zeros(3,lattice.totF);
            lattice.faces.nodes = zeros(nodesperface,lattice.totF);
            lattice.faces.edges = zeros(edgesperface,lattice.totF);
            %%----------------------------------------------------- lattice edges ------------------------------------------------------------------------------------------------
            for j=1:Ndim2
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + 1;
              
                    edgeindex1 = k + (j-1)*(Ndim1-1);
      
                    lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                    lattice.edges.nodes(2,edgeindex1) = nodeindex2;
      
                    lattice.nodes.edges(1,nodeindex1) = edgeindex1;
                    lattice.nodes.edges(3,nodeindex2) = edgeindex1;
      
                    lattice.edges.id(edgeindex1) = edgeindex1;
                    lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,index1) + lattice.nodes.meshcoordinates(1,index2))/2;
                    lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,index1) + lattice.nodes.meshcoordinates(2,index2))/2;
                    lattice.edges.centroid(3,edgeindex1) = 0;
                end
            end
            for j=1:Ndim2-1
                for k=1:Ndim1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
      
                    edgeindex1 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
       
                    lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                    lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                    lattice.nodes.edges(2,nodeindex1) = edgeindex1;
                    lattice.nodes.edges(4,nodeindex2) = edgeindex1;
       
                    lattice.edges.id(edgeindex1) = edgeindex1;
                    lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,index1) + lattice.nodes.meshcoordinates(1,index2))/2;
                    lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,index1) + lattice.nodes.meshcoordinates(2,index2))/2;
                    lattice.edges.centroid(3,edgeindex1) = 0;
                end
            end
            %%----------------------------------------------------- lattice faces ------------------------------------------------------------------------------------------------
            for j=1:Ndim2-1
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
                    nodeindex3 = k + (j-1)*Ndim1 + Ndim1 + 1;
                    nodeindex4 = k + (j-1)*Ndim1 + 1;

                    edgeindex1 = k + (j-1)*(Ndim1-1);
                    edgeindex2 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
                    edgeindex3 = k + (j-1)*(Ndim1-1) + (Ndim1-1);
                    edgeindex4 = k + (j-1)*(Ndim1-1) + 1;
       
                    faceindex1 = k + (j-1)*(Ndim1-1);
       
                    lattice.faces.nodes(1,faceindex1) = index1;
                    lattice.faces.nodes(2,faceindex1) = index2;
                    lattice.faces.nodes(3,faceindex1) = index3;
                    lattice.faces.nodes(4,faceindex1) = index4;
       
                    lattice.nodes.faces(1,nodeindex1) = faceindex1;
                    lattice.nodes.faces(4,nodeindex2) = faceindex1;
                    lattice.nodes.faces(3,nodeindex3) = faceindex1;
                    lattice.nodes.faces(2,nodeindex4) = faceindex1;
       
                    lattice.faces.id(faceindex1) = faceindex1;
                    lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                    lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                    lattice.faces.centroid(3,faceindex1) = 0; 
                    
                    lattice.faces.edges(1,faceindex1) = edgeindex1;
                    lattice.faces.edges(2,faceindex1) = edgeindex2;
                    lattice.faces.edges(3,faceindex1) = edgeindex3;
                    lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                    lattice.edges.faces(1,edgeindex1) = faceindex1;
                    lattice.edges.faces(1,edgeindex2) = faceindex1;
                    lattice.edges.faces(2,edgeindex3) = faceindex1;
                    lattice.edges.faces(2,edgeindex4) = faceindex1; 
                end
            end
        elseif edgeflag
            lattice.edgeflag = 1;
            lattice.totE = (lattice.Ndim1-1)*lattice.Ndim2 + (lattice.Ndim2-1)*lattice.Ndim1;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
                disp(' ')
            end
            lattice.nodesperedge = 2;
            lattice.edgespernode = 4;
            lattice.nodes.edges = zeros(edgespernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(nodesperedge,lattice.totE);
            %%----------------------------------------------------- lattice edges ------------------------------------------------------------------------------------------------
            for j=1:Ndim2
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + 1;
              
                    edgeindex1 = k + (j-1)*(Ndim1-1);
      
                    lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                    lattice.edges.nodes(2,edgeindex1) = nodeindex2;
      
                    lattice.nodes.edges(1,nodeindex1) = edgeindex1;
                    lattice.nodes.edges(3,nodeindex2) = edgeindex1;
      
                    lattice.edges.id(edgeindex1) = edgeindex1;
                    lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,index1) + lattice.nodes.meshcoordinates(1,index2))/2;
                    lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,index1) + lattice.nodes.meshcoordinates(2,index2))/2;
                    lattice.edges.centroid(3,edgeindex1) = 0;
                end
            end
            for j=1:Ndim2-1
                for k=1:Ndim1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
      
                    edgeindex1 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
       
                    lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                    lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                    lattice.nodes.edges(2,nodeindex1) = edgeindex1;
                    lattice.nodes.edges(4,nodeindex2) = edgeindex1;
       
                    lattice.edges.id(edgeindex1) = edgeindex1;
                    lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,index1) + lattice.nodes.meshcoordinates(1,index2))/2;
                    lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,index1) + lattice.nodes.meshcoordinates(2,index2))/2;
                    lattice.edges.centroid(3,edgeindex1) = 0;
                end
            end
        end
%-----------------------------------------------------------------------------3D--------------------------------------------------------------------------------------------
    case 3
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('lattice GENERATION')
            disp(' ')
            disp('Meshing of a cuboid with brick elements - 3D')
            disp(' ')
            fprintf('Lower boundary along the first direction x0 = %d\n',lattice.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',lattice.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',lattice.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',lattice.dim2max)
            fprintf('Lower boundary along the second direction z0 = %d\n',lattice.dim3min)
            fprintf('Upper boundary along the second direction zN = %6.2f\n',lattice.dim3max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',lattice.totN)
        end
        lattice.nodes.id = zeros(1,lattice.totN);
        lattice.nodes.meshcoordinates = zeros(3,lattice.totN);  
        lattice.deltadim1 = deltadim1;
        lattice.deltadim2 = deltadim2;
        lattice.deltadim3 = deltadim3;
        for i=1:lattice.Ndim3
            for j=1:lattice.Ndim2
                for k=1:lattice.Ndim1
                    lattice.nodes.meshcoordinates(1,k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2) = lattice.dim1min + lattice.deltadim1*(k-1);
                    lattice.nodes.meshcoordinates(2,k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2) = lattice.dim2min + lattice.deltadim2*(j-1); 
                    lattice.nodes.meshcoordinates(3,k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2) = lattice.dim3min + lattice.deltadim3*(i-1);
                    lattice.nodes.id(k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2) = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2;
                end
            end
        end
        lattice.coordlines.dim1free = zeros(lattice.Nlinesdim1,3*lattice.Ndim2*lattice.Ndim3);
        lattice.coordlines.dim2free = zeros(lattice.Nlinesdim2,3*lattice.Ndim1*lattice.Ndim3);
        lattice.coordlines.dim3free = zeros(lattice.Nlinesdim3,3*lattice.Ndim1*lattice.Ndim2);
        for i=1:lattice.Ndim3
            for j=1:lattice.Ndim2
                coord2 = lattice.dim2min + lattice.deltadim2*(j-1);
                coord3 = lattice.dim3min + lattice.deltadim3*(i-1);
                deltalinedim1 = (lattice.dim1max - lattice.dim1min)/(lattice.Nlinesdim1 - 1);
                for l=1:lattice.Nlinesdim1
                    lattice.coordlines.dim1free(l,3*(j-1)+3*lattice.Ndim2*(i-1) + 1) = lattice.dim1max + deltalinedim1*(l-1);
                    lattice.coordlines.dim1free(l,3*(j-1)+3*lattice.Ndim2*(i-1) + 2) = coord2; 
                    lattice.coordlines.dim1free(l,3*(j-1)+3*lattice.Ndim2*(i-1) + 3) = coord3;
                end
            end
        end
        for i=1:lattice.Ndim3
            for k=1:lattice.Ndim1
                coord1 = lattice.dim1min + lattice.deltadim1*(k-1);
                coord3 = lattice.dim3min + lattice.deltadim3*(i-1);
                deltalinedim2 = (lattice.dim2max - lattice.dim2min)/(lattice.Nlinesdim2 - 1);
                for l=1:lattice.Nlinesdim2
                    lattice.coordlines.dim2free(l,3*(k-1)+3*lattice.Ndim1*(i-1) + 1) = coord1;
                    lattice.coordlines.dim2free(l,3*(k-1)+3*lattice.Ndim1*(i-1) + 2) = lattice.dim2max + deltalinedim2*(l-1); 
                    lattice.coordlines.dim2free(l,3*(k-1)+3*lattice.Ndim1*(i-1) + 3) = coord3;
                end
            end
        end
        for j=1:lattice.Ndim2
            for k=1:lattice.Ndim1
                coord1 = lattice.dim1min + lattice.deltadim1*(k-1);
                coord2 = lattice.dim2min + lattice.deltadim2*(j-1);
                deltalinedim3 = (lattice.dim3max - lattice.dim3min)/(lattice.Nlinesdim3 - 1);
                for l=1:lattice.Nlinesdim3
                    lattice.coordlines.dim3free(l,3*(k-1)+3*lattice.Ndim1*(j-1) + 1) = coord1;
                    lattice.coordlines.dim3free(l,3*(k-1)+3*lattice.Ndim1*(j-1) + 2) = coord2; 
                    lattice.coordlines.dim3free(l,3*(k-1)+3*lattice.Ndim1*(j-1) + 3) = lattice.dim3max + deltalinedim3*(l-1);
                end
            end
        end
        if (edgeflag && faceflag && cellflag)
            lattice.edgeflag = 1;
            lattice.faceflag = 1;
            lattice.cellflag = 1;
            lattice.totE = (lattice.Ndim1-1)*lattice.Ndim2*lattice.Ndim3 + (lattice.Ndim2-1)*lattice.Ndim1*lattice.Ndim3 + (lattice.Ndim3-1)*lattice.Ndim1*lattice.Ndim2;
            lattice.totF = (lattice.Ndim1-1)*(lattice.Ndim2-1)*lattice.Ndim3 + (lattice.Ndim2-1)*(lattice.Ndim3-1)*lattice.Ndim1 + (lattice.Ndim1-1)*(lattice.Ndim3-1)*lattice.Ndim2;
            lattice.totC = (lattice.Ndim1-1)*(lattice.Ndim2-1)*(lattice.Ndim3-1);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
                fprintf('Total number of faces F = %6.2f\n',lattice.totF)
                fprintf('Total number of cells C = %6.2f\n',lattice.totC)
                disp(' ')
            end
            lattice.nodesperedge = 2;
            lattice.nodesperface = 4;
            lattice.nodespercell = 8;
            lattice.edgespernode = 6;
            lattice.edgesperface = 4;
            lattice.edgespercell = 12;
            lattice.facespernode = 12;
            lattice.facesperedge = 4;
            lattice.facespercell = 6;
            lattice.cellspernode = 8;
            lattice.cellsperedge = 4;
            lattice.cellsperface = 2;
            lattice.nodes.edges = zeros(lattice.edgespernode,lattice.totN);
            lattice.nodes.faces = zeros(lattice.facespernode,lattice.totN);
            lattice.nodes.cells = zeros(lattice.cellspernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(lattice.nodesperedge,lattice.totE);
            lattice.edges.faces = zeros(lattice.facesperedge,lattice.totE);
            lattice.edges.cells = zeros(lattice.cellsperedge,lattice.totE);
            lattice.faces.id = zeros(1,lattice.totF);
            lattice.faces.centroid = zeros(3,lattice.totF);
            lattice.faces.nodes = zeros(lattice.nodesperface,lattice.totF);
            lattice.faces.edges = zeros(lattice.edgesperface,lattice.totF);
            lattice.faces.cells = zeros(lattice.cellsperface,lattice.totF);
            lattice.cells.id = zeros(1,lattice.totC);
            lattice.cells.centroid = zeros(3,lattice.totC);
            lattice.cells.nodes = zeros(lattice.nodespercell,lattice.totC);
            lattice.cells.edges = zeros(lattice.edgespercell,lattice.totC);
            lattice.cells.faces = zeros(lattice.facespercell,lattice.totC);
            %-------------------------------------->edges
            for i=1:lattice.Ndim3
                for j=1:lattice.Ndim2
                    for k=1:lattice.Ndim1-1
                        nodeindex1 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2;
                        nodeindex2 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(lattice.Ndim1-1) + (i-1)*(lattice.Ndim1-1)*lattice.Ndim2;
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(2,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(1,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(4,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(3,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(6,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(5,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            %-------------------------------------->faces
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1 + 1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1);
                        edgeindex4 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + 1;
       
                        faceindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(1,nodeindex1) = faceindex1;
                        lattice.nodes.faces(4,nodeindex2) = faceindex1;
                        lattice.nodes.faces(3,nodeindex3) = faceindex1;
                        lattice.nodes.faces(2,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(1,edgeindex1) = faceindex1;
                        lattice.edges.faces(1,edgeindex2) = faceindex1;
                        lattice.edges.faces(2,edgeindex3) = faceindex1;
                        lattice.edges.faces(2,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + 1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1)*Ndim2;
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
       
                        faceindex1 = (Ndim1-1)*(Ndim2-1)*Ndim3 + k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                                                
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(5,nodeindex1) = faceindex1;
                        lattice.nodes.faces(8,nodeindex2) = faceindex1;
                        lattice.nodes.faces(7,nodeindex3) = faceindex1;
                        lattice.nodes.faces(6,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(3,edgeindex1) = faceindex1;
                        lattice.edges.faces(1,edgeindex2) = faceindex1;
                        lattice.edges.faces(4,edgeindex3) = faceindex1;
                        lattice.edges.faces(2,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + Ndim1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;

                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex3 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + Ndim1*(Ndim2-1);
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
       
                        faceindex1 = (Ndim1-1)*(Ndim2-1)*Ndim3 + (Ndim1-1)*(Ndim3-1)*Ndim2 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(9,nodeindex1) = faceindex1;
                        lattice.nodes.faces(12,nodeindex2) = faceindex1;
                        lattice.nodes.faces(11,nodeindex3) = faceindex1;
                        lattice.nodes.faces(10,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(3,edgeindex1) = faceindex1;
                        lattice.edges.faces(3,edgeindex2) = faceindex1;
                        lattice.edges.faces(4,edgeindex3) = faceindex1;
                        lattice.edges.faces(4,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            %-------------------------------------->cells
            for i=1:Ndim3-1
                for j=1:Ndim2-1
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1 + 1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        nodeindex5 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex6 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + Ndim1;
                        nodeindex7 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + Ndim1 + 1;
                        nodeindex8 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + 1;
       
                        edgeindex1 =  k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 =  (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex3 =  k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1);
                        edgeindex4 =  (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + 1;
                        edgeindex5 =  (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex6 =  (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex7 =  (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1 + 1;
                        edgeindex8 =  (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        edgeindex9 =  k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1)*Ndim2;
                        edgeindex10 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + Ndim1*(Ndim2-1);
                        edgeindex11 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1) + (Ndim1-1)*Ndim2;
                        edgeindex12 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + 1 + Ndim1*(Ndim2-1);
       
                        faceindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        faceindex2 = (Ndim1-1)*(Ndim2-1)*Ndim3 + k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        faceindex3 = (Ndim1-1)*(Ndim2-1)*Ndim3 + (Ndim1-1)*(Ndim3-1)*Ndim2 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        faceindex4 = (Ndim1-1)*(Ndim2-1)*Ndim3 + k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1);
                        faceindex5 = (Ndim1-1)*(Ndim2-1)*Ndim3 + (Ndim1-1)*(Ndim3-1)*Ndim2 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + 1;
                        faceindex6 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1) + (Ndim1-1)*(Ndim2-1);
       
                        cellindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        
                        lattice.cells.nodes(1,cellindex1) = nodeindex1;
                        lattice.cells.nodes(2,cellindex1) = nodeindex2;
                        lattice.cells.nodes(3,cellindex1) = nodeindex3;
                        lattice.cells.nodes(4,cellindex1) = nodeindex4;
                        lattice.cells.nodes(5,cellindex1) = nodeindex5;
                        lattice.cells.nodes(6,cellindex1) = nodeindex6;
                        lattice.cells.nodes(7,cellindex1) = nodeindex7;
                        lattice.cells.nodes(8,cellindex1) = nodeindex8;
       
                        lattice.nodes.cells(1,nodeindex1) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(4,nodeindex2) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(3,nodeindex3) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(2,nodeindex4) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(5,nodeindex5) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(8,nodeindex6) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(7,nodeindex7) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
                        lattice.nodes.cells(6,nodeindex8) = k + j*(Ndim1-1) + i*(Ndim1-1)*(Ndim2-1);
       
                        lattice.cells.id(cellindex1) = cellindex1;
                        lattice.cells.centroid(1,cellindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4) + lattice.nodes.meshcoordinates(1,nodeindex5) + lattice.nodes.meshcoordinates(1,nodeindex6) + lattice.nodes.meshcoordinates(1,nodeindex7) + lattice.nodes.meshcoordinates(1,nodeindex8))/8;
                        lattice.cells.centroid(2,cellindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4) + lattice.nodes.meshcoordinates(2,nodeindex5) + lattice.nodes.meshcoordinates(2,nodeindex6) + lattice.nodes.meshcoordinates(2,nodeindex7) + lattice.nodes.meshcoordinates(2,nodeindex8))/8;
                        lattice.cells.centroid(3,cellindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4) + lattice.nodes.meshcoordinates(3,nodeindex5) + lattice.nodes.meshcoordinates(3,nodeindex6) + lattice.nodes.meshcoordinates(3,nodeindex7) + lattice.nodes.meshcoordinates(3,nodeindex8))/8;
       
                        lattice.cells.edges(1,cellindex1) = edgeindex1;
                        lattice.cells.edges(2,cellindex1) = edgeindex2;
                        lattice.cells.edges(3,cellindex1) = edgeindex3;
                        lattice.cells.edges(4,cellindex1) = edgeindex4;
                        lattice.cells.edges(5,cellindex1) = edgeindex5;
                        lattice.cells.edges(6,cellindex1) = edgeindex6;
                        lattice.cells.edges(7,cellindex1) = edgeindex7;
                        lattice.cells.edges(8,cellindex1) = edgeindex8;
                        lattice.cells.edges(9,cellindex1) = edgeindex9;
                        lattice.cells.edges(10,cellindex1) = edgeindex10;
                        lattice.cells.edges(11,cellindex1) = edgeindex11;
                        lattice.cells.edges(12,cellindex1) = edgeindex12;
       
                        lattice.edges.cells(1,edgeindex1) = cellindex1;
                        lattice.edges.cells(1,edgeindex2) = cellindex1;
                        lattice.edges.cells(2,edgeindex3) = cellindex1;
                        lattice.edges.cells(2,edgeindex4) = cellindex1;
                        lattice.edges.cells(1,edgeindex5) = cellindex1;
                        lattice.edges.cells(4,edgeindex6) = cellindex1;
                        lattice.edges.cells(3,edgeindex7) = cellindex1;
                        lattice.edges.cells(2,edgeindex8) = cellindex1;
                        lattice.edges.cells(4,edgeindex9) = cellindex1;
                        lattice.edges.cells(4,edgeindex10) = cellindex1;
                        lattice.edges.cells(3,edgeindex11) = cellindex1;
                        lattice.edges.cells(3,edgeindex12) = cellindex1;
       
                        lattice.cells.faces(1,cellindex1) = faceindex1;
                        lattice.cells.faces(2,cellindex1) = faceindex2;
                        lattice.cells.faces(3,cellindex1) = faceindex3;
                        lattice.cells.faces(4,cellindex1) = faceindex4;
                        lattice.cells.faces(5,cellindex1) = faceindex5;
                        lattice.cells.faces(6,cellindex1) = faceindex6;
       
                        lattice.faces.cells(1,faceindex1) = cellindex1;
                        lattice.faces.cells(1,faceindex2) = cellindex1;
                        lattice.faces.cells(1,faceindex3) = cellindex1;
                        lattice.faces.cells(2,faceindex4) = cellindex1;
                        lattice.faces.cells(2,faceindex5) = cellindex1;
                        lattice.faces.cells(2,faceindex6) = cellindex1;
                    end
                end
            end
        elseif (edgeflag && faceflag)
            lattice.edgeflag = 1;
            lattice.faceflag = 1;
            lattice.totE = (lattice.Ndim1-1)*lattice.Ndim2*lattice.Ndim3 + (lattice.Ndim2-1)*lattice.Ndim1*lattice.Ndim3 + (lattice.Ndim3-1)*lattice.Ndim1*lattice.Ndim2;
            lattice.totF = (lattice.Ndim1-1)*(lattice.Ndim2-1)*lattice.Ndim3 + (lattice.Ndim2-1)*(lattice.Ndim3-1)*lattice.Ndim1 + (lattice.Ndim1-1)*(lattice.Ndim3-1)*lattice.Ndim2;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
                fprintf('Total number of faces F = %6.2f\n',lattice.totF)
                disp(' ')
            end
            lattice.nodesperedge = 2;
            lattice.nodesperface = 4;
            lattice.edgespernode = 6;
            lattice.edgesperface = 4;
            lattice.facespernode = 12;
            lattice.facesperedge = 4;
            lattice.nodes.edges = zeros(lattice.edgespernode,lattice.totN);
            lattice.nodes.faces = zeros(lattice.facespernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(lattice.nodesperedge,lattice.totE);
            lattice.edges.faces = zeros(lattice.facesperedge,lattice.totE);
            lattice.faces.id = zeros(1,lattice.totF);
            lattice.faces.centroid = zeros(3,lattice.totF);
            lattice.faces.nodes = zeros(lattice.nodesperface,lattice.totF);
            lattice.faces.edges = zeros(lattice.edgesperface,lattice.totF);
            %-------------------------------------->edges
            for i=1:Ndim3
                for j=1:Ndim2
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(2,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(1,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(4,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(3,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(6,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(5,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            %-------------------------------------->faces
            for i=1:lattice.Ndim3
                for j=1:lattice.Ndim2
                    for k=1:lattice.Ndim1-1
                        nodeindex1 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2;
                        nodeindex2 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2 + lattice.Ndim1;
                        nodeindex3 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2 + lattice.Ndim1 + 1;
                        nodeindex4 = k + (j-1)*lattice.Ndim1 + (i-1)*lattice.Ndim1*lattice.Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1);
                        edgeindex4 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + 1;
       
                        faceindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(1,nodeindex1) = faceindex1;
                        lattice.nodes.faces(4,nodeindex2) = faceindex1;
                        lattice.nodes.faces(3,nodeindex3) = faceindex1;
                        lattice.nodes.faces(2,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(1,edgeindex1) = faceindex1;
                        lattice.edges.faces(1,edgeindex2) = faceindex1;
                        lattice.edges.faces(2,edgeindex3) = faceindex1;
                        lattice.edges.faces(2,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            for i=1:lattice.Ndim3
                for j=1:lattice.Ndim2-1
                    for k=1:lattice.Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + 1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1)*Ndim2;
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
       
                        faceindex1 = (Ndim1-1)*(Ndim2-1)*Ndim3 + k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(5,nodeindex1) = faceindex1;
                        lattice.nodes.faces(8,nodeindex2) = faceindex1;
                        lattice.nodes.faces(7,nodeindex3) = faceindex1;
                        lattice.nodes.faces(6,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(3,edgeindex1) = faceindex1;
                        lattice.edges.faces(1,edgeindex2) = faceindex1;
                        lattice.edges.faces(4,edgeindex3) = faceindex1;
                        lattice.edges.faces(2,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + Ndim1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;

                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex3 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + Ndim1*(Ndim2-1);
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
       
                        faceindex1 = (Ndim1-1)*(Ndim2-1)*Ndim3 + (Ndim1-1)*(Ndim3-1)*Ndim2 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                                                
                        lattice.faces.nodes(1,faceindex1) = nodeindex1;
                        lattice.faces.nodes(2,faceindex1) = nodeindex2;
                        lattice.faces.nodes(3,faceindex1) = nodeindex3;
                        lattice.faces.nodes(4,faceindex1) = nodeindex4;
       
                        lattice.nodes.faces(9,nodeindex1) = faceindex1;
                        lattice.nodes.faces(12,nodeindex2) = faceindex1;
                        lattice.nodes.faces(11,nodeindex3) = faceindex1;
                        lattice.nodes.faces(10,nodeindex4) = faceindex1;
       
                        lattice.faces.id(faceindex1) = faceindex1;
                        lattice.faces.centroid(1,faceindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1) + lattice.nodes.meshcoordinates(1,nodeindex2) + lattice.nodes.meshcoordinates(1,nodeindex3) + lattice.nodes.meshcoordinates(1,nodeindex4))/4;
                        lattice.faces.centroid(2,faceindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1) + lattice.nodes.meshcoordinates(2,nodeindex2) + lattice.nodes.meshcoordinates(2,nodeindex3) + lattice.nodes.meshcoordinates(2,nodeindex4))/4;
                        lattice.faces.centroid(3,faceindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1) + lattice.nodes.meshcoordinates(3,nodeindex2) + lattice.nodes.meshcoordinates(3,nodeindex3) + lattice.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        lattice.faces.edges(1,faceindex1) = edgeindex1;
                        lattice.faces.edges(2,faceindex1) = edgeindex2;
                        lattice.faces.edges(3,faceindex1) = edgeindex3;
                        lattice.faces.edges(4,faceindex1) = edgeindex4;
       
                        lattice.edges.faces(3,edgeindex1) = faceindex1;
                        lattice.edges.faces(3,edgeindex2) = faceindex1;
                        lattice.edges.faces(4,edgeindex3) = faceindex1;
                        lattice.edges.faces(4,edgeindex4) = faceindex1;
                   
                    end
                end
            end
        elseif edgeflag
            lattice.edgeflag = 1;
            lattice.totE = (lattice.Ndim1-1)*lattice.Ndim2*lattice.Ndim3 + (lattice.Ndim2-1)*lattice.Ndim1*lattice.Ndim3 + (lattice.Ndim3-1)*lattice.Ndim1*lattice.Ndim2;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',lattice.totE)
                disp(' ')
            end
            lattice.nodesperedge = 2;
            lattice.edgespernode = 6;
            lattice.nodes.edges = zeros(lattice.edgespernode,lattice.totN);
            lattice.edges.id = zeros(1,lattice.totE);
            lattice.edges.centroid = zeros(3,lattice.totE);
            lattice.edges.nodes = zeros(lattice.nodesperedge,lattice.totE);
            %-------------------------------------->edges
            for i=1:Ndim3
                for j=1:Ndim2
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(2,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(1,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(4,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(3,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        lattice.edges.nodes(1,edgeindex1) = nodeindex1;
                        lattice.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        lattice.nodes.edges(6,nodeindex1) = edgeindex1;
                        lattice.nodes.edges(5,nodeindex1) = edgeindex1;
       
                        lattice.edges.id(edgeindex1) = edgeindex1;
                        lattice.edges.centroid(1,edgeindex1) = (lattice.nodes.meshcoordinates(1,nodeindex1)+lattice.nodes.meshcoordinates(1,nodeindex2))/2;
                        lattice.edges.centroid(2,edgeindex1) = (lattice.nodes.meshcoordinates(2,nodeindex1)+lattice.nodes.meshcoordinates(2,nodeindex2))/2;
                        lattice.edges.centroid(3,edgeindex1) = (lattice.nodes.meshcoordinates(3,nodeindex1)+lattice.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return
  
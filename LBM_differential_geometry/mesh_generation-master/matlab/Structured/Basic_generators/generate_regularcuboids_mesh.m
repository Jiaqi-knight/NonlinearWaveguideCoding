function[mesh]=generate_regularcuboids_mesh(D,dim1min,dim1max,Ndim1,dim2min,dim2max,Ndim2,dim3min,dim3max,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 27th, 2014
%    Last update: April 8th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

mesh.D = D;
mesh.dim1min = dim1min;
mesh.dim1max = dim1max;
mesh.Ndim1 = Ndim1;
mesh.dim2min = dim2min;
mesh.dim2max = dim2max;
mesh.Ndim2 = Ndim2;
mesh.dim3min = dim3min;
mesh.dim3max = dim3max;
mesh.Ndim3 = Ndim3;
mesh.Nlinesdim1 = Nlinesdim1;
mesh.Nlinesdim2 = Nlinesdim2;
mesh.Nlinesdim3 = Nlinesdim3;
mesh.totN = Ndim1*Ndim2*Ndim3;
mesh.totE = 0;
mesh.totF = 0;
mesh.totC = 0;
mesh.edgeflag = edgeflag;
mesh.faceflag = faceflag;
mesh.cellflag = cellflag;
   
switch D
%-----------------------------------------------------------------------------1D--------------------------------------------------------------------------------------------
    case 1
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('MESH GENERATION')
            disp(' ')
            disp('Generation of equally spaced intervals along a line - 1D')
            fprintf('Start point x0 = %d\n',mesh.dim1min)
            fprintf('End point xN = %6.2f\n',mesh.dim1max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);  
        mesh.deltadim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Ndim1 - 1);
        for i=1:Ndim1
            mesh.nodes.meshcoordinates(1,i) = mesh.dim1min + mesh.deltadim1*(i-1);
            mesh.nodes.meshcoordinates(2,i) = 0; 
            mesh.nodes.meshcoordinates(3,i) = 0;
            mesh.nodes.id(i) = i;
        end
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3);
        coord2 = 0;
        coord3 = 0;
        deltalinedim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Nlinesdim1 - 1);
        for l=1:mesh.Nlinesdim1
            mesh.coordlines.dim1free(l,1) = mesh.dim1min + deltalinedim1*(l-1);
            mesh.coordlines.dim1free(l,2) = coord2; 
            mesh.coordlines.dim1free(l,3) = coord3;
        end
        if edgeflag || faceflag || cellflag
            mesh.edgeflag = 1;
            mesh.edgespernode = 2;
            mesh.nodesperedge = 2;
            mesh.totE = (mesh.Ndim1-1);
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
            end
            for i=1:mesh.Ndim1-1
               nodeindex1 = i;
               nodeindex2 = i + 1;
               edgeindex1 = i;
               mesh.nodes.edges(2,nodeindex1) = edgeindex1;
               mesh.nodes.edges(1,nodeindex2) = edgeindex1;
               mesh.edges.id(i) = edgeindex1;
               mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
               mesh.edges.centroid(2,i) = 0;
               mesh.edges.centroid(3,i) = 0;
               mesh.edges.nodes(1,i) = nodeindex1;
               mesh.edges.nodes(2,i) = nodeindex2;
            end
        end
%-----------------------------------------------------------------------------2D--------------------------------------------------------------------------------------------
    case 2
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('MESH GENERATION')
            disp(' ')
            disp('Meshing of a rectangular domain with quadrilateral elements - 2D')
            fprintf('Lower boundary along the first direction x0 = %d\n',mesh.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',mesh.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',mesh.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',mesh.dim2max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);  
        mesh.deltadim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Ndim1 - 1);
        mesh.deltadim2 = (mesh.dim2max - mesh.dim2min)/(mesh.Ndim2 - 1);
        for i=1:Ndim2
            for j=1:Ndim1
                mesh.nodes.meshcoordinates(1,j + (i-1)*Ndim1) = mesh.dim1min + mesh.deltadim1*(j-1);
                mesh.nodes.meshcoordinates(2,j + (i-1)*Ndim1) = mesh.dim2min + mesh.deltadim2*(i-1); 
                mesh.nodes.meshcoordinates(3,j + (i-1)*Ndim1) = 0;
                mesh.nodes.id(j + i*Ndim1) = j + (i-1)*Ndim1;
            end
        end
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3*mesh.Ndim2);
        mesh.coordlines.dim2free = zeros(mesh.Nlinesdim2,3*mesh.Ndim1);
        for j=1:mesh.Ndim2
            coord2 = mesh.dim2min + mesh.deltadim2*(j-1);
            coord3 = 0;
            deltalinedim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Nlinesdim1 - 1);
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,3*(j-1) + 1) = mesh.dim1min + deltalinedim1*(l-1);
                mesh.coordlines.dim1free(l,3*(j-1) + 2) = coord2; 
                mesh.coordlines.dim1free(l,3*(j-1) + 3) = coord3;
            end
        end
        for k=1:mesh.Ndim1
            coord1 = mesh.dim1min + mesh.deltadim1*(k-1);
            coord3 = 0;
            deltalinedim2 = (mesh.dim2max - mesh.dim2min)/(mesh.Nlinesdim2 - 1);
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,3*(k-1) + 1) = coord1;
                mesh.coordlines.dim2free(l,3*(k-1) + 2) = mesh.dim2min + deltalinedim2*(l-1); 
                mesh.coordlines.dim2free(l,3*(k-1) + 3) = coord3;
            end
        end
        if (edgeflag && faceflag) || (edgeflag && cellflag)
            mesh.edgeflag = 1;
            mesh.faceflag = 1;
            mesh.totE = (mesh.Ndim1-1)*mesh.Ndim2 + (mesh.Ndim2-1)*mesh.Ndim1;
            mesh.totF = (mesh.Ndim1-1)*(mesh.Ndim2-1);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
                fprintf('Total number of faces F = %6.2f\n',mesh.totF)
                disp(' ')
            end
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.edgespernode = 4;
            mesh.edgesperface = 4;
            mesh.facespernode = 4;
            mesh.facesperedge = 2;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            %%----------------------------------------------------- mesh edges ------------------------------------------------------------------------------------------------
            for j=1:Ndim2
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + 1;
              
                    edgeindex1 = k + (j-1)*(Ndim1-1);
      
                    mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                    mesh.edges.nodes(2,edgeindex1) = nodeindex2;
      
                    mesh.nodes.edges(1,nodeindex1) = edgeindex1;
                    mesh.nodes.edges(3,nodeindex2) = edgeindex1;
      
                    mesh.edges.id(edgeindex1) = edgeindex1;
                    mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                    mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                    mesh.edges.centroid(3,edgeindex1) = 0;
                end
            end
            for j=1:Ndim2-1
                for k=1:Ndim1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
      
                    edgeindex1 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
       
                    mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                    mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                    mesh.nodes.edges(2,nodeindex1) = edgeindex1;
                    mesh.nodes.edges(4,nodeindex2) = edgeindex1;
       
                    mesh.edges.id(edgeindex1) = edgeindex1;
                    mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                    mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                    mesh.edges.centroid(3,edgeindex1) = 0;
                end
            end
            %%----------------------------------------------------- mesh faces ------------------------------------------------------------------------------------------------
            for j=1:Ndim2-1
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
                    nodeindex3 = k + (j-1)*Ndim1 + Ndim1 + 1;
                    nodeindex4 = k + (j-1)*Ndim1 + 1;

                    edgeindex1 = k + (j-1)*(Ndim1-1);
                    edgeindex2 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
                    edgeindex3 = k + (j-1)*(Ndim1-1) + (Ndim1-1);
                    edgeindex4 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1 + 1;
       
                    faceindex1 = k + (j-1)*(Ndim1-1);
       
                    mesh.faces.nodes(1,faceindex1) = nodeindex1;
                    mesh.faces.nodes(2,faceindex1) = nodeindex4;
                    mesh.faces.nodes(3,faceindex1) = nodeindex3;
                    mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                    mesh.nodes.faces(1,nodeindex1) = faceindex1;
                    mesh.nodes.faces(4,nodeindex2) = faceindex1;
                    mesh.nodes.faces(3,nodeindex3) = faceindex1;
                    mesh.nodes.faces(2,nodeindex4) = faceindex1;
       
                    mesh.faces.id(faceindex1) = faceindex1;
                    mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                    mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                    mesh.faces.centroid(3,faceindex1) = 0; 
                    
                    mesh.faces.edges(1,faceindex1) = edgeindex1;
                    mesh.faces.edges(2,faceindex1) = edgeindex2;
                    mesh.faces.edges(3,faceindex1) = edgeindex3;
                    mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                    mesh.edges.faces(1,edgeindex1) = faceindex1;
                    mesh.edges.faces(1,edgeindex2) = faceindex1;
                    mesh.edges.faces(2,edgeindex3) = faceindex1;
                    mesh.edges.faces(2,edgeindex4) = faceindex1; 
                end
            end
        elseif edgeflag
            mesh.edgeflag = 1;
            mesh.totE = (mesh.Ndim1-1)*mesh.Ndim2 + (mesh.Ndim2-1)*mesh.Ndim1;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
                disp(' ')
            end
            mesh.nodesperedge = 2;
            mesh.edgespernode = 4;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            %%----------------------------------------------------- mesh edges ------------------------------------------------------------------------------------------------
            for j=1:Ndim2
                for k=1:Ndim1-1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + 1;
              
                    edgeindex1 = k + (j-1)*(Ndim1-1);
      
                    mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                    mesh.edges.nodes(2,edgeindex1) = nodeindex2;
      
                    mesh.nodes.edges(1,nodeindex1) = edgeindex1;
                    mesh.nodes.edges(3,nodeindex2) = edgeindex1;
      
                    mesh.edges.id(edgeindex1) = edgeindex1;
                    mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                    mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                    mesh.edges.centroid(3,edgeindex1) = 0;
                end
            end
            for j=1:Ndim2-1
                for k=1:Ndim1
                    nodeindex1 = k + (j-1)*Ndim1;
                    nodeindex2 = k + (j-1)*Ndim1 + Ndim1;
      
                    edgeindex1 = (Ndim1-1)*Ndim2 + k + (j-1)*Ndim1;
       
                    mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                    mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                    mesh.nodes.edges(2,nodeindex1) = edgeindex1;
                    mesh.nodes.edges(4,nodeindex2) = edgeindex1;
       
                    mesh.edges.id(edgeindex1) = edgeindex1;
                    mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                    mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                    mesh.edges.centroid(3,edgeindex1) = 0;
                end
            end
        end
%-----------------------------------------------------------------------------3D--------------------------------------------------------------------------------------------
    case 3
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('MESH GENERATION')
            disp(' ')
            disp('Meshing of a cuboid with brick elements - 3D')
            disp(' ')
            fprintf('Lower boundary along the first direction x0 = %d\n',mesh.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',mesh.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',mesh.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',mesh.dim2max)
            fprintf('Lower boundary along the second direction z0 = %d\n',mesh.dim3min)
            fprintf('Upper boundary along the second direction zN = %6.2f\n',mesh.dim3max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);  
        mesh.deltadim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Ndim1 - 1);
        mesh.deltadim2 = (mesh.dim2max - mesh.dim2min)/(mesh.Ndim2 - 1);
        mesh.deltadim3 = (mesh.dim3max - mesh.dim3min)/(mesh.Ndim3 - 1);
        for i=1:mesh.Ndim3
            for j=1:mesh.Ndim2
                for k=1:mesh.Ndim1
                    mesh.nodes.meshcoordinates(1,k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2) = mesh.dim1min + mesh.deltadim1*(k-1);
                    mesh.nodes.meshcoordinates(2,k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2) = mesh.dim2min + mesh.deltadim2*(j-1); 
                    mesh.nodes.meshcoordinates(3,k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2) = mesh.dim3min + mesh.deltadim3*(i-1);
                    mesh.nodes.id(k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2) = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2;
                end
            end
        end
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3*mesh.Ndim2*mesh.Ndim3);
        mesh.coordlines.dim2free = zeros(mesh.Nlinesdim2,3*mesh.Ndim1*mesh.Ndim3);
        mesh.coordlines.dim3free = zeros(mesh.Nlinesdim3,3*mesh.Ndim1*mesh.Ndim2);
        for i=1:mesh.Ndim3
            for j=1:mesh.Ndim2
                coord2 = mesh.dim2min + mesh.deltadim2*(j-1);
                coord3 = mesh.dim3min + mesh.deltadim3*(i-1);
                deltalinedim1 = (mesh.dim1max - mesh.dim1min)/(mesh.Nlinesdim1 - 1);
                for l=1:mesh.Nlinesdim1
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1) = mesh.dim1min + deltalinedim1*(l-1);
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2) = coord2; 
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3) = coord3;
                end
            end
        end
        for i=1:mesh.Ndim3
            for k=1:mesh.Ndim1
                coord1 = mesh.dim1min + mesh.deltadim1*(k-1);
                coord3 = mesh.dim3min + mesh.deltadim3*(i-1);
                deltalinedim2 = (mesh.dim2max - mesh.dim2min)/(mesh.Nlinesdim2 - 1);
                for l=1:mesh.Nlinesdim2
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1) = coord1;
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2) = mesh.dim2min + deltalinedim2*(l-1); 
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3) = coord3;
                end
            end
        end
        for j=1:mesh.Ndim2
            for k=1:mesh.Ndim1
                coord1 = mesh.dim1min + mesh.deltadim1*(k-1);
                coord2 = mesh.dim2min + mesh.deltadim2*(j-1);
                deltalinedim3 = (mesh.dim3max - mesh.dim3min)/(mesh.Nlinesdim3 - 1);
                for l=1:mesh.Nlinesdim3
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1) = coord1;
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2) = coord2; 
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3) = mesh.dim3min + deltalinedim3*(l-1);
                end
            end
        end
        if (edgeflag && faceflag && cellflag)
            mesh.edgeflag = 1;
            mesh.faceflag = 1;
            mesh.cellflag = 1;
            mesh.totE = (mesh.Ndim1-1)*mesh.Ndim2*mesh.Ndim3 + (mesh.Ndim2-1)*mesh.Ndim1*mesh.Ndim3 + (mesh.Ndim3-1)*mesh.Ndim1*mesh.Ndim2;
            mesh.totF = (mesh.Ndim1-1)*(mesh.Ndim2-1)*mesh.Ndim3 + (mesh.Ndim2-1)*(mesh.Ndim3-1)*mesh.Ndim1 + (mesh.Ndim1-1)*(mesh.Ndim3-1)*mesh.Ndim2;
            mesh.totC = (mesh.Ndim1-1)*(mesh.Ndim2-1)*(mesh.Ndim3-1);
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
                fprintf('Total number of faces F = %6.2f\n',mesh.totF)
                fprintf('Total number of cells C = %6.2f\n',mesh.totC)
                disp(' ')
            end
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.nodespercell = 8;
            mesh.edgespernode = 6;
            mesh.edgesperface = 4;
            mesh.edgespercell = 12;
            mesh.facespernode = 12;
            mesh.facesperedge = 4;
            mesh.facespercell = 6;
            mesh.cellspernode = 8;
            mesh.cellsperedge = 4;
            mesh.cellsperface = 2;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.nodes.cells = zeros(mesh.cellspernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.edges.cells = zeros(mesh.cellsperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            mesh.faces.cells = zeros(mesh.cellsperface,mesh.totF);
            mesh.cells.id = zeros(1,mesh.totC);
            mesh.cells.centroid = zeros(3,mesh.totC);
            mesh.cells.nodes = zeros(mesh.nodespercell,mesh.totC);
            mesh.cells.edges = zeros(mesh.edgespercell,mesh.totC);
            mesh.cells.faces = zeros(mesh.facespercell,mesh.totC);
            %-------------------------------------->edges
            for i=1:mesh.Ndim3
                for j=1:mesh.Ndim2
                    for k=1:mesh.Ndim1-1
                        nodeindex1 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2;
                        nodeindex2 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(mesh.Ndim1-1) + (i-1)*(mesh.Ndim1-1)*mesh.Ndim2;
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(2,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(1,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(4,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(3,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(6,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(5,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
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
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + 1;
       
                        faceindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(1,nodeindex1) = faceindex1;
                        mesh.nodes.faces(4,nodeindex2) = faceindex1;
                        mesh.nodes.faces(3,nodeindex3) = faceindex1;
                        mesh.nodes.faces(2,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(1,edgeindex1) = faceindex1;
                        mesh.edges.faces(1,edgeindex2) = faceindex1;
                        mesh.edges.faces(2,edgeindex3) = faceindex1;
                        mesh.edges.faces(2,edgeindex4) = faceindex1;
                   
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
                                                
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(5,nodeindex1) = faceindex1;
                        mesh.nodes.faces(8,nodeindex2) = faceindex1;
                        mesh.nodes.faces(7,nodeindex3) = faceindex1;
                        mesh.nodes.faces(6,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(3,edgeindex1) = faceindex1;
                        mesh.edges.faces(1,edgeindex2) = faceindex1;
                        mesh.edges.faces(4,edgeindex3) = faceindex1;
                        mesh.edges.faces(2,edgeindex4) = faceindex1;
                   
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
                        
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(9,nodeindex1) = faceindex1;
                        mesh.nodes.faces(12,nodeindex2) = faceindex1;
                        mesh.nodes.faces(11,nodeindex3) = faceindex1;
                        mesh.nodes.faces(10,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(3,edgeindex1) = faceindex1;
                        mesh.edges.faces(3,edgeindex2) = faceindex1;
                        mesh.edges.faces(4,edgeindex3) = faceindex1;
                        mesh.edges.faces(4,edgeindex4) = faceindex1;
                   
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
                        
                        mesh.cells.nodes(1,cellindex1) = nodeindex1;
                        mesh.cells.nodes(2,cellindex1) = nodeindex4;
                        mesh.cells.nodes(3,cellindex1) = nodeindex3;
                        mesh.cells.nodes(4,cellindex1) = nodeindex2;
                        mesh.cells.nodes(5,cellindex1) = nodeindex5;
                        mesh.cells.nodes(6,cellindex1) = nodeindex8;
                        mesh.cells.nodes(7,cellindex1) = nodeindex7;
                        mesh.cells.nodes(8,cellindex1) = nodeindex6;
       
                        mesh.nodes.cells(1,nodeindex1) = cellindex1;
                        mesh.nodes.cells(4,nodeindex2) = cellindex1;
                        mesh.nodes.cells(3,nodeindex3) = cellindex1;
                        mesh.nodes.cells(2,nodeindex4) = cellindex1;
                        mesh.nodes.cells(5,nodeindex5) = cellindex1;
                        mesh.nodes.cells(8,nodeindex6) = cellindex1;
                        mesh.nodes.cells(7,nodeindex7) = cellindex1;
                        mesh.nodes.cells(6,nodeindex8) = cellindex1;
       
                        mesh.cells.id(cellindex1) = cellindex1;
                        mesh.cells.centroid(1,cellindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4) + mesh.nodes.meshcoordinates(1,nodeindex5) + mesh.nodes.meshcoordinates(1,nodeindex6) + mesh.nodes.meshcoordinates(1,nodeindex7) + mesh.nodes.meshcoordinates(1,nodeindex8))/8;
                        mesh.cells.centroid(2,cellindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4) + mesh.nodes.meshcoordinates(2,nodeindex5) + mesh.nodes.meshcoordinates(2,nodeindex6) + mesh.nodes.meshcoordinates(2,nodeindex7) + mesh.nodes.meshcoordinates(2,nodeindex8))/8;
                        mesh.cells.centroid(3,cellindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4) + mesh.nodes.meshcoordinates(3,nodeindex5) + mesh.nodes.meshcoordinates(3,nodeindex6) + mesh.nodes.meshcoordinates(3,nodeindex7) + mesh.nodes.meshcoordinates(3,nodeindex8))/8;
       
                        mesh.cells.edges(1,cellindex1) = edgeindex1;
                        mesh.cells.edges(2,cellindex1) = edgeindex2;
                        mesh.cells.edges(3,cellindex1) = edgeindex3;
                        mesh.cells.edges(4,cellindex1) = edgeindex4;
                        mesh.cells.edges(5,cellindex1) = edgeindex5;
                        mesh.cells.edges(6,cellindex1) = edgeindex6;
                        mesh.cells.edges(7,cellindex1) = edgeindex7;
                        mesh.cells.edges(8,cellindex1) = edgeindex8;
                        mesh.cells.edges(9,cellindex1) = edgeindex9;
                        mesh.cells.edges(10,cellindex1) = edgeindex10;
                        mesh.cells.edges(11,cellindex1) = edgeindex11;
                        mesh.cells.edges(12,cellindex1) = edgeindex12;
       
                        mesh.edges.cells(1,edgeindex1) = cellindex1;
                        mesh.edges.cells(1,edgeindex2) = cellindex1;
                        mesh.edges.cells(2,edgeindex3) = cellindex1;
                        mesh.edges.cells(2,edgeindex4) = cellindex1;
                        mesh.edges.cells(1,edgeindex5) = cellindex1;
                        mesh.edges.cells(4,edgeindex6) = cellindex1;
                        mesh.edges.cells(3,edgeindex7) = cellindex1;
                        mesh.edges.cells(2,edgeindex8) = cellindex1;
                        mesh.edges.cells(4,edgeindex9) = cellindex1;
                        mesh.edges.cells(4,edgeindex10) = cellindex1;
                        mesh.edges.cells(3,edgeindex11) = cellindex1;
                        mesh.edges.cells(3,edgeindex12) = cellindex1;
       
                        mesh.cells.faces(1,cellindex1) = faceindex1;
                        mesh.cells.faces(2,cellindex1) = faceindex2;
                        mesh.cells.faces(3,cellindex1) = faceindex3;
                        mesh.cells.faces(4,cellindex1) = faceindex4;
                        mesh.cells.faces(5,cellindex1) = faceindex5;
                        mesh.cells.faces(6,cellindex1) = faceindex6;
       
                        mesh.faces.cells(1,faceindex1) = cellindex1;
                        mesh.faces.cells(1,faceindex2) = cellindex1;
                        mesh.faces.cells(1,faceindex3) = cellindex1;
                        mesh.faces.cells(2,faceindex4) = cellindex1;
                        mesh.faces.cells(2,faceindex5) = cellindex1;
                        mesh.faces.cells(2,faceindex6) = cellindex1;
                    end
                end
            end
        elseif (edgeflag && faceflag)
            mesh.edgeflag = 1;
            mesh.faceflag = 1;
            mesh.totE = (mesh.Ndim1-1)*mesh.Ndim2*mesh.Ndim3 + (mesh.Ndim2-1)*mesh.Ndim1*mesh.Ndim3 + (mesh.Ndim3-1)*mesh.Ndim1*mesh.Ndim2;
            mesh.totF = (mesh.Ndim1-1)*(mesh.Ndim2-1)*mesh.Ndim3 + (mesh.Ndim2-1)*(mesh.Ndim3-1)*mesh.Ndim1 + (mesh.Ndim1-1)*(mesh.Ndim3-1)*mesh.Ndim2;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
                fprintf('Total number of faces F = %6.2f\n',mesh.totF)
                disp(' ')
            end
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.edgespernode = 6;
            mesh.edgesperface = 4;
            mesh.facespernode = 12;
            mesh.facesperedge = 4;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            %-------------------------------------->edges
            for i=1:Ndim3
                for j=1:Ndim2
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(2,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(1,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(4,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(3,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(6,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(5,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            %-------------------------------------->faces
            for i=1:mesh.Ndim3
                for j=1:mesh.Ndim2
                    for k=1:mesh.Ndim1-1
                        nodeindex1 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2;
                        nodeindex2 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2 + mesh.Ndim1;
                        nodeindex3 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2 + mesh.Ndim1 + 1;
                        nodeindex4 = k + (j-1)*mesh.Ndim1 + (i-1)*mesh.Ndim1*mesh.Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1);
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1) + 1;
       
                        faceindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*(Ndim2-1);
                        
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(1,nodeindex1) = faceindex1;
                        mesh.nodes.faces(4,nodeindex2) = faceindex1;
                        mesh.nodes.faces(3,nodeindex3) = faceindex1;
                        mesh.nodes.faces(2,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(1,edgeindex1) = faceindex1;
                        mesh.edges.faces(1,edgeindex2) = faceindex1;
                        mesh.edges.faces(2,edgeindex3) = faceindex1;
                        mesh.edges.faces(2,edgeindex4) = faceindex1;
                   
                    end
                end
            end
            for i=1:mesh.Ndim3
                for j=1:mesh.Ndim2-1
                    for k=1:mesh.Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        nodeindex3 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2 + 1;
                        nodeindex4 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;

                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        edgeindex2 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        edgeindex3 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2 + (Ndim1-1)*Ndim2;
                        edgeindex4 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
       
                        faceindex1 = (Ndim1-1)*(Ndim2-1)*Ndim3 + k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
                        
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(5,nodeindex1) = faceindex1;
                        mesh.nodes.faces(8,nodeindex2) = faceindex1;
                        mesh.nodes.faces(7,nodeindex3) = faceindex1;
                        mesh.nodes.faces(6,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(3,edgeindex1) = faceindex1;
                        mesh.edges.faces(1,edgeindex2) = faceindex1;
                        mesh.edges.faces(4,edgeindex3) = faceindex1;
                        mesh.edges.faces(2,edgeindex4) = faceindex1;
                   
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
                                                
                        mesh.faces.nodes(1,faceindex1) = nodeindex1;
                        mesh.faces.nodes(2,faceindex1) = nodeindex4;
                        mesh.faces.nodes(3,faceindex1) = nodeindex3;
                        mesh.faces.nodes(4,faceindex1) = nodeindex2;
       
                        mesh.nodes.faces(9,nodeindex1) = faceindex1;
                        mesh.nodes.faces(12,nodeindex2) = faceindex1;
                        mesh.nodes.faces(11,nodeindex3) = faceindex1;
                        mesh.nodes.faces(10,nodeindex4) = faceindex1;
       
                        mesh.faces.id(faceindex1) = faceindex1;
                        mesh.faces.centroid(1,faceindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1) + mesh.nodes.meshcoordinates(1,nodeindex2) + mesh.nodes.meshcoordinates(1,nodeindex3) + mesh.nodes.meshcoordinates(1,nodeindex4))/4;
                        mesh.faces.centroid(2,faceindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1) + mesh.nodes.meshcoordinates(2,nodeindex2) + mesh.nodes.meshcoordinates(2,nodeindex3) + mesh.nodes.meshcoordinates(2,nodeindex4))/4;
                        mesh.faces.centroid(3,faceindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1) + mesh.nodes.meshcoordinates(3,nodeindex2) + mesh.nodes.meshcoordinates(3,nodeindex3) + mesh.nodes.meshcoordinates(3,nodeindex4))/4;
       
                        mesh.faces.edges(1,faceindex1) = edgeindex1;
                        mesh.faces.edges(2,faceindex1) = edgeindex2;
                        mesh.faces.edges(3,faceindex1) = edgeindex3;
                        mesh.faces.edges(4,faceindex1) = edgeindex4;
       
                        mesh.edges.faces(3,edgeindex1) = faceindex1;
                        mesh.edges.faces(3,edgeindex2) = faceindex1;
                        mesh.edges.faces(4,edgeindex3) = faceindex1;
                        mesh.edges.faces(4,edgeindex4) = faceindex1;
                   
                    end
                end
            end
        elseif edgeflag
            mesh.edgeflag = 1;
            mesh.totE = (mesh.Ndim1-1)*mesh.Ndim2*mesh.Ndim3 + (mesh.Ndim2-1)*mesh.Ndim1*mesh.Ndim3 + (mesh.Ndim3-1)*mesh.Ndim1*mesh.Ndim2;
            if printflag
                fprintf('Total number of edges E = %6.2f\n',mesh.totE)
                disp(' ')
            end
            mesh.nodesperedge = 2;
            mesh.edgespernode = 6;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            %-------------------------------------->edges
            for i=1:Ndim3
                for j=1:Ndim2
                    for k=1:Ndim1-1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + 1;
                        edgeindex1 = k + (j-1)*(Ndim1-1) + (i-1)*(Ndim1-1)*Ndim2;
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(2,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(1,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
            for i=1:Ndim3
                for j=1:Ndim2-1
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*(Ndim2-1);
       
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(4,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(3,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;                   
                    end
                end
            end
            for i=1:Ndim3-1
                for j=1:Ndim2
                    for k=1:Ndim1
                        nodeindex1 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
                        nodeindex2 = k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2 + Ndim1*Ndim2;
                        edgeindex1 = (Ndim1-1)*Ndim2*Ndim3+(Ndim2-1)*Ndim1*Ndim3 + k + (j-1)*Ndim1 + (i-1)*Ndim1*Ndim2;
      
                        mesh.edges.nodes(1,edgeindex1) = nodeindex1;
                        mesh.edges.nodes(2,edgeindex1) = nodeindex2;
       
                        mesh.nodes.edges(6,nodeindex1) = edgeindex1;
                        mesh.nodes.edges(5,nodeindex2) = edgeindex1;
       
                        mesh.edges.id(edgeindex1) = edgeindex1;
                        mesh.edges.centroid(1,edgeindex1) = (mesh.nodes.meshcoordinates(1,nodeindex1)+mesh.nodes.meshcoordinates(1,nodeindex2))/2;
                        mesh.edges.centroid(2,edgeindex1) = (mesh.nodes.meshcoordinates(2,nodeindex1)+mesh.nodes.meshcoordinates(2,nodeindex2))/2;
                        mesh.edges.centroid(3,edgeindex1) = (mesh.nodes.meshcoordinates(3,nodeindex1)+mesh.nodes.meshcoordinates(3,nodeindex2))/2;
                    end
                end
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return
  